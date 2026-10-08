"""Depletion and excess: when could the sea be seen? -- step 24 companion.

Verifies ``docs/analysis/sea_depletion.md``.  Step 23's Proposition Q1 says
that E never reads the sea: the sea enters an event only as its source or
its sink.  So the sea can be seen only through a rule that reads it.  This
demo asks which rules there are, what they cost the observable, and how deep
a sea makes them invisible.  Everything runs on step 16's mesh, with its
compensated kernel, generalised to any potential.

Two ledgers.  ``SeaLedger`` is step 16's ``Ledger`` for any potential: u+
and u- are transported separately and spectrally.  ``MinimalLedger``
keeps the ensemble minimal (garbage collection, Proposition Q8), so
u+ = E+ and u- = E- at every step and only E is transported: the ledger
then carries no transport artefact (Part F says why that matters).

A  Proposition V1, admissible negativity is free.  For a potential of
   degree two or less the compensated residual kernel vanishes identically,
   so free and harmonic evolution make no events and spend no sea pair,
   whatever the state.  A cat state in a harmonic trap: no events, the sea
   untouched.  The same cat in an anharmonic well: events, and more of them
   than a single packet, because its fringes are bodies too.
B  Proposition V2, the sea's debt is the negativity the dynamics makes:
   in the minimal ledger D = dM- - dN_tr/2 exactly, dN_tr being the change
   of ||E||_1 that the discrete transport itself makes.
C  Proposition V3, the floor.  A sea of depth beta B from which an emission
   must draw a pair at its parent's cell, never below zero.  E is unchanged
   bit for bit for beta >= lambda* (the unfloored run's worst deficit) and
   deviates below, through a spurious force, an energy drift and a
   shortfall of negativity.  Absorptive-first and emissive realisations,
   under (S) and (S').
D  Proposition V4, granular supply.  At the physical density the event
   aperture dp * 2 y_max holds one sea pair on average (Q9(c)), so an
   integer sea is empty there with probability exp(-n).  The expected E
   under that rule against the depth beta.
E  lambda* over long runs (open item Q-SP2) and against the reach.
R  (not in the default) the transport's share of dM- in Part B against
   the grid: Part B's minimal ledger on finer meshes, about 25 minutes.
G  Addendum: pair collisions, the amended (S′) (Propositions V8-V10).  A
   W-null motion of the sea's pair density on top of (S′) streaming:
   Definition (C), Focus/Defocus between aligned pairs (mass action with
   detailed balance), against a reversible Volterra sea and a symmetric
   pair-diffusion control with the kernel's channels and rates.  SymPy
   identities (Volterra conserves F, collisions decrease H), lambda*(t),
   a rate scan, other cases and the reach.  Sub-parts with --g-subs:
   I exact and invariants, M main comparison, S rate scan, R other cases,
   X reach ladder.
F  Demo defect: step 16's ledger books transport ringing as sea traffic.
   In a harmonic trap, with no events at all, it inflates the body count
   and debits the sea; the share of the published Eckart run's traffic
   that is ringing.

Run as::

    WPMW_OUTPUT=... PYTHONPATH=src python3 -u src/demo_sea_depletion.py \\
        [--parts ABCDEFG] [--quick] [--t-long 96] [--g-subs IMSRX]

About forty minutes for all parts on one core (``--parts`` lets them run
in parallel); ``--quick`` halves the run lengths, for testing only.
Parts C, D, E and G write CSV files through ``output_path``, and the floor,
long-run and collision figures are drawn from those files, so parts may run in separate
processes; run any part last to redraw.
"""

from __future__ import annotations

import argparse
import csv
import time

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import demo_sea_population_equilibrium as d16
from wpmwlib.wpmw_utils import output_path, docs_path

HBAR = d16.HBAR
B = d16.B                       # 1/(pi hbar), the crystal shift 2/h


def banner(title: str) -> None:
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72, flush=True)


def save_fig(fig, name: str) -> None:
    fig.savefig(output_path(name), dpi=150, bbox_inches="tight")
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=150, bbox_inches="tight")
    plt.close(fig)


def write_csv(name: str, header, rows) -> None:
    with open(output_path(name), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def rel_l2(a, b):
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


# ----------------------------------------------------------------------
# Potentials
# ----------------------------------------------------------------------
def eckart(v0=1.0, a=1.0):
    """The step 15 barrier, V0 sech^2(r/a)."""
    return lambda r: v0 / np.cosh(r / a) ** 2


def well(v0=4.0, a=2.0):
    """An attractive Poeschl-Teller well, -V0 sech^2(r/a): bounded and
    anharmonic, so a displaced packet dephases and makes negativity."""
    return lambda r: -v0 / np.cosh(r / a) ** 2


def harmonic(omega=1.0):
    return lambda r: 0.5 * omega ** 2 * r ** 2


def free():
    return lambda r: 0.0 * r


# ----------------------------------------------------------------------
class SeaLedger(d16.Ledger):
    """Step 16's ledger for any potential, with supply rules for the sea.

    The symbols are built exactly as in ``d16.Ledger`` and ``d16.Horizon``:
    raised-cosine window first, on both symbols, then the compensation by
    the kernel's own discrete first moment (spec 3.2 and 4.1).  ``y_h=None``
    puts the window at the grid's Nyquist reach, as ``Ledger`` does, and the
    symbols are then identical to ``Ledger``'s for the Eckart barrier.
    """

    def __init__(self, Vf, n_r=128, r_half=20.0, n_p=64, dp=0.25, y_h=None):
        self.Vf = Vf
        self.v0, self.a, self.n_p, self.dp = 0.0, 1.0, n_p, dp
        self.r = -r_half + 2.0 * r_half * np.arange(n_r) / n_r
        self.dr = self.r[1] - self.r[0]
        self.kr = 2.0 * np.pi * np.fft.fftfreq(n_r, d=self.dr)
        self.p = dp * (np.fft.fftfreq(n_p, d=1.0 / n_p).astype(int) + 0.5)
        s = 2.0 * np.pi * np.fft.fftfreq(n_p, d=dp)
        self.y = HBAR * s / 2.0
        self.y_max = float(np.abs(self.y).max())
        self.area = self.dr * self.dp
        self.xi = np.round(np.fft.fftfreq(n_p, d=1.0) * n_p) * dp
        nyq = n_p // 2
        rr, yy = self.r[:, None], self.y[None, :]
        m_full = (1j / HBAR) * (Vf(rr + yy) - Vf(rr - yy))
        s_sym = np.broadcast_to(s, m_full.shape).copy()
        m_full[:, nyq] = 0.0
        s_sym[:, nyq] = 0.0

        def first_moment(sym):
            return np.real((self.xi * np.fft.ifft(sym, axis=-1)).sum(axis=-1))

        if y_h is None:
            w = np.cos(np.pi * self.y / (2.0 * self.y_max)) ** 2
            self.y_h = self.y_max
        else:
            w = np.where(np.abs(self.y) <= y_h,
                         np.cos(np.pi * self.y / (2.0 * y_h)) ** 2, 0.0)
            self.y_h = y_h
        m_full, s_sym = m_full * w[None, :], s_sym * w[None, :]
        self.dv_eff = first_moment(m_full) / first_moment(1j * s_sym)
        m_res = m_full - 1j * self.dv_eff[:, None] * s_sym
        self.k = np.real(np.fft.ifft(m_res, axis=1))
        self.sym_e = np.fft.fft(self.k, axis=1)
        self.sym_a = np.real(np.fft.fft(np.abs(self.k), axis=1))
        self.gamma_tot = self.sym_a[:, 0].copy()
        self.cls = 1j * self.dv_eff[:, None] * s[None, :]
        self.Vr = Vf(self.r)

    # -- events with a supply rule -------------------------------------
    def channels_supply(self, up, um, S, dt, rule="none", nu=1.0, est=None,
                        absorb=True, ring=None, smin=None):
        """``d16.Ledger.channels`` with a supply rule for emission.

        Each channel pair (q, -q) is ONE event; absorptive if a partner of
        the right species sits at BOTH daughters (and ``absorb``), else
        emissive, splitting a pair at the parent's cell.
        rule 'none'    : the sea is a diagnostic ledger, as in step 16.
        rule 'floor'   : an emission cannot draw the parent's cell below 0.
        rule 'poisson' : the parent's aperture holds Poisson(nu S/B) pairs;
                         an emission is realised with probability
                         1 - exp(-nu S/B).  Mean-field expectation.
        A blocked emission is lost: its two depositions do not happen.
        ``est`` accumulates, for an unfloored run, the emission demand and
        the share an integer sea would have blocked at each nu.  ``ring``
        accumulates the negative part of the absorptive allocation, which
        only a ringing (negative) partner field can produce (Part F).
        ``smin`` tracks the floor's margin: the lowest (S - Em)/B over
        cells where an emission is due, S read before the channel's
        credits, which is exactly what rule 'floor' compares.  At depth
        beta B the floor binds iff beta - 1 + margin_1 < 0 (Q1: the sea at
        depth beta is the sea at depth 1 plus a constant), so
        lambda* = 1 - margin_1 is the depth at which it starts to bind.
        """
        n_abs = n_emi = blocked = 0.0
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
                if absorb:
                    capA = np.where(t > 0, np.roll(um, -q, axis=1),
                                    np.roll(up, -q, axis=1))
                    capB = np.where(t > 0, np.roll(up, q, axis=1),
                                    np.roll(um, q, axis=1))
                    A = np.minimum(D, np.minimum(capA, capB))
                    if ring is not None:
                        ring[0] += float(np.maximum(-A, 0.0).sum())
                else:
                    A = np.zeros_like(D)
                Em = D - A
                if smin is not None:
                    due = Em > 0.0
                    if due.any():
                        smin[0] = min(smin[0],
                                      float((S - Em)[due].min()) / B)
                if est is not None:
                    est["demand"] += float(Em.sum())
                    s_ = np.maximum(S, 0.0) / B
                    for j, v in enumerate(est["nus"]):
                        est["blocked"][j] += float((Em * np.exp(-v * s_)).sum())
                if rule == "floor":
                    Er = np.minimum(Em, np.maximum(S, 0.0))
                elif rule == "poisson":
                    Er = Em * (1.0 - np.exp(-nu * np.maximum(S, 0.0) / B))
                else:
                    Er = Em
                blocked += float((Em - Er).sum())
                Em = Er
                n_abs += float(A.sum())
                n_emi += float(Em.sum())
                aq, am = np.roll(A, q, axis=1), np.roll(A, -q, axis=1)
                um -= np.where(t > 0, aq, 0.0)
                up -= np.where(t > 0, 0.0, aq)
                up -= np.where(t > 0, am, 0.0)
                um -= np.where(t > 0, 0.0, am)
                S += A
                eq, em_ = np.roll(Em, q, axis=1), np.roll(Em, -q, axis=1)
                up += np.where(t > 0, eq, 0.0)
                um += np.where(t > 0, 0.0, eq)
                um += np.where(t > 0, em_, 0.0)
                up += np.where(t > 0, 0.0, em_)
                S -= Em
                np.maximum(up, 0.0, out=up)
                np.maximum(um, 0.0, out=um)
        return up, um, S, n_abs, n_emi, blocked

    # -- step 16's two-field run, for Part F ----------------------------
    def run16(self, e0, t_max, dt, sea_force=True):
        """``d16.Ledger.run_events`` (absorptive) for this potential, with
        the ringing counter.  Bit for bit the published ledger."""
        self.sea_force = sea_force
        up = np.maximum(e0, 0.0).copy()
        um = np.maximum(-e0, 0.0).copy()
        S = np.full_like(e0, B)
        ring = [0.0]
        tot = np.zeros(2)
        mn = 1.0
        for _ in range(int(round(t_max / dt))):
            up, um, S = (self.stream(up, .5 * dt), self.stream(um, .5 * dt),
                         self.stream_sea(S, .5 * dt))
            up, um, S, a, e_, _b = self.channels_supply(up, um, S, dt,
                                                        ring=ring)
            tot += (a, e_)
            up, um, S = (self.stream(up, .5 * dt), self.stream(um, .5 * dt),
                         self.stream_sea(S, .5 * dt))
            mn = min(mn, float((S / B).min()))
        return dict(e=up - um, N=float((up + um).sum()) * self.area,
                    S=S, min_s=mn, ring=ring[0] * self.area,
                    n_abs=tot[0] * self.area, n_emi=tot[1] * self.area)

    # -- diagnostics -------------------------------------------------------
    def neg_mass(self, e):
        return float(np.maximum(-e, 0.0).sum()) * self.area

    def l1(self, e):
        return float(np.abs(e).sum()) * self.area

    def moments(self, e):
        a = self.area
        p = self.p[None, :]
        return dict(norm=float(e.sum()) * a,
                    p=float((e * p).sum()) * a,
                    H=float((e * (0.5 * p ** 2 + self.Vr[:, None])).sum()) * a,
                    Mneg=self.neg_mass(e),
                    T=float(e[self.r > 0.0].sum()) * a)


# ----------------------------------------------------------------------
class MinimalLedger:
    """The ledger with an immediate contact sink (Proposition Q8).

    With garbage collection every cell holds one species only after each
    step, so u+ = E+ and u- = E-.  Exact transport carries disjoint fields
    to disjoint fields, so it keeps them minimal; only E need be moved, and
    E is the signed field the spectral transport is accurate for.  Each
    step: half stream E and S; events from u+ = E+, u- = E- with
    ``SeaLedger.channels_supply``; contact recombination c = min(u+, u-),
    credited to S at the cell; E = u+ - u-; half stream.  The discrete
    transport does not conserve ||E||_1 exactly; that change, dN_tr, is
    booked in its own column and never charged to the sea.
    """

    def __init__(self, mesh: SeaLedger):
        self.m = mesh
        self.area = mesh.area

    def run(self, e0, t_max, dt, sea_force=True, beta=1.0, rule="none",
            nu=1.0, absorb=True, every=None, nus=None, ref=True):
        m = self.m
        m.sea_force = sea_force
        a = self.area
        E = e0.copy()
        S = np.full_like(e0, beta * B)
        S0 = float(S.sum()) * a
        M0 = m.neg_mass(e0)
        n_tr = 0.0
        r_e = e0.copy()
        tot = np.zeros(4)                       # abs, emi, blocked, contact
        deficit, excess = 0.0, 0.0
        smin = [beta]                           # the floor's margin, in B
        est = None if nus is None else dict(
            nus=list(nus), demand=0.0, blocked=np.zeros(len(nus)))
        tr = []
        n_steps = int(round(t_max / dt))
        k_every = None if every is None else max(1, int(round(every / dt)))
        for step in range(n_steps + 1):
            if k_every and (step % k_every == 0 or step == n_steps):
                tr.append((step * dt, S0 - float(S.sum()) * a,
                           m.neg_mass(E) - M0, 0.5 * n_tr,
                           beta - float(S.min()) / B,
                           float(S.max()) / B - beta, tot[2] * a,
                           beta - smin[0]))
            if step == n_steps:
                break
            for half in (0, 1):
                if half == 1:
                    up, um = np.maximum(E, 0.0), np.maximum(-E, 0.0)
                    up, um, S, x1, x2, x3 = m.channels_supply(
                        up, um, S, dt, rule, nu, est, absorb, smin=smin)
                    c = np.minimum(up, um)
                    S += c
                    E = up - um
                    tot += (x1, x2, x3, float(c.sum()))
                    deficit = beta - smin[0]
                    excess = max(excess, float(S.max()) / B - beta)
                n0 = m.l1(E)
                E = m.stream(E, 0.5 * dt)
                n_tr += m.l1(E) - n0
                S = m.stream_sea(S, 0.5 * dt)
            if ref:
                r_e = m.qle_step(r_e, dt)
        n_ev = tot[0] + tot[1]
        out = dict(e=E, S=S, ref=r_e, deficit=deficit, excess=excess,
                   f=tot[0] / max(n_ev, 1e-30), n_ev=n_ev * a,
                   blocked=tot[2] * a, contact=tot[3] * a,
                   blocked_frac=tot[2] / max(tot[1] + tot[2], 1e-30),
                   N=m.l1(E), n_tr=n_tr, trace=np.array(tr), est=est)
        if ref:
            out["fid"] = rel_l2(E, r_e)
        return out


# ----------------------------------------------------------------------
# States
# ----------------------------------------------------------------------
def coherent(run, x0, p0, sigma):
    """Gaussian Wigner function of a minimum-uncertainty packet."""
    sp = HBAR / (2.0 * sigma)
    w = np.exp(-((run.r[:, None] - x0) ** 2) / (2.0 * sigma ** 2)
               - ((run.p[None, :] - p0) ** 2) / (2.0 * sp ** 2))
    return w / (w.sum() * run.area)


def cat(run, x0, sigma):
    """Wigner function of the even cat |x0> + |-x0>, Gaussians of width
    sigma: two lobes and the fringe term 2 G(x, p) cos(2 x0 p / hbar)."""
    r, p = run.r[:, None], run.p[None, :]

    def g(x):
        return np.exp(-x ** 2 / (2.0 * sigma ** 2)
                      - 2.0 * sigma ** 2 * p ** 2 / HBAR ** 2)

    w = g(r - x0) + g(r + x0) + 2.0 * g(r) * np.cos(2.0 * x0 * p / HBAR)
    return w / (w.sum() * run.area)


def summit_packet(run):
    """The step 16 benchmark: Eckart summit, p = 1, sigma_r = 1."""
    return d16.gaussian(run, 0.0, 1.0, 1.0)


def well_packet(run):
    """A packet displaced in the well, at rest: it oscillates and dephases."""
    return coherent(run, 2.0, 0.0, 0.6)


CASES = (("Eckart", eckart, summit_packet, 8.0),
         ("PT well", well, well_packet, 12.0))
MODES = (("absorb-first", True, True), ("absorb-first", False, True),
         ("emissive", False, False))        # (name, sea_force, absorb)


def mode_label(force, absorb):
    return f"{'(S)' if force else '(S′)'} {'abs' if absorb else 'emi'}"


# ----------------------------------------------------------------------
# A.  Proposition V1
# ----------------------------------------------------------------------
def part_a(q):
    banner("A  Proposition V1: admissible negativity is free")
    print("  The residual kernel is what is left of the Moyal symbol")
    print("  (i/hbar)[V(x+y) - V(x-y)] once the force 2 y V'_eff(x) is")
    print("  removed.  For deg V <= 2 nothing is left.\n")
    print(f"  {'potential':>10} {'max|V_eff_prime|':>17} {'max|K|':>11}"
          f" {'ratio':>11}")
    for name, Vf in (("free", free()), ("harmonic", harmonic()),
                     ("Eckart", eckart()), ("PT well", well())):
        run = SeaLedger(Vf)
        kmax = float(np.abs(run.k[:, 1:]).max())
        dvm = float(np.abs(run.dv_eff).max())
        rat = f"{kmax / dvm:11.3e}" if dvm > 0 else f"{'--':>11}"
        print(f"  {name:>10} {dvm:17.4e} {kmax:11.3e} {rat}")

    T = np.pi if q else 2.0 * np.pi
    dt = 0.02
    out = {}
    print(f"\n  Minimal ledger, (S'), x0 = 3, sigma = 1/sqrt(2), T = {T:.3f},"
          f" dt = {dt}:")
    print(f"  {'case':>16} {'M-(0)':>8} {'M-(T)':>8} {'events':>10}"
          f" {'per t':>8} {'debt D':>9} {'lambda*':>9} {'|E-QLE|':>9}")
    rows = (("cat, harmonic", harmonic(), "cat"),
            ("cat, PT well", well(), "cat"),
            ("packet, PT well", well(), "packet"))
    for name, Vf, kind in rows:
        run = SeaLedger(Vf)
        e0 = (cat(run, 3.0, 1.0 / np.sqrt(2.0)) if kind == "cat"
              else coherent(run, 3.0, 0.0, 1.0 / np.sqrt(2.0)))
        o = MinimalLedger(run).run(e0, T, dt, sea_force=False)
        out[name] = (run, e0, o)
        dbt = float((B - o["S"]).sum()) * run.area
        print(f"  {name:>16} {run.neg_mass(e0):8.4f}"
              f" {run.neg_mass(o['e']):8.4f} {o['n_ev']:10.3e}"
              f" {o['n_ev'] / T:8.3f} {dbt:9.5f} {o['deficit']:9.4f}"
              f" {o['fid']:9.2e}")
    print("\n  In the trap the cat's negative mass is carried by negatons")
    print("  from the start and rotates rigidly; no event happens and no")
    print("  sea pair is touched.  In the well the same cat makes events at")
    print("  a higher rate than a single packet of the same width and")
    print("  offset: its fringes are bodies, and every body is a parent.")
    return out


# ----------------------------------------------------------------------
# B.  Proposition V2
# ----------------------------------------------------------------------
def part_b(q):
    banner("B  Proposition V2: the sea's debt is the negativity made")
    print("  Minimal ledger.  Theorem S7 event by event and the contact")
    print("  sink give dS = -dN/2 for everything but transport; N = ||E||_1")
    print("  and the norm is conserved, so")
    print("      D = S0 - S_tot = dM- - dN_tr/2   exactly.\n")
    out = {}
    print(f"  {'case':>8} {'motion':>9} {'t':>5} {'D':>9} {'dM-':>9}"
          f" {'dN_tr/2':>9} {'residual':>9}")
    for name, Vf, prep, T in CASES:
        if q:
            T *= 0.5
        run = SeaLedger(Vf())
        e0 = prep(run)
        for force, absorb in ((True, True), (False, True)):
            o = MinimalLedger(run).run(e0, T, 0.02, sea_force=force,
                                       absorb=absorb, every=T / 4.0,
                                       ref=False)
            out[(name, force)] = o
            for row in o["trace"][1:]:
                t, D, dm, ntr2 = row[:4]
                print(f"  {name[:8]:>8} {mode_label(force, absorb):>9}"
                      f" {t:5.1f} {D:9.5f} {dm:9.5f} {ntr2:+9.5f}"
                      f" {D - dm + ntr2:9.1e}")
    print("\n  The residual is round-off.  The sea pays exactly for the")
    print("  negative mass the dynamics makes, less what the discrete")
    print("  transport makes or destroys on its own.")
    return out


# ----------------------------------------------------------------------
# C.  Proposition V3: the floor
# ----------------------------------------------------------------------
def part_c(q):
    banner("C  Proposition V3: the floor (minimal ledger)")
    print("  Sea depth beta B.  An emission draws its pair at the parent's")
    print("  cell and cannot draw below zero; a blocked emission is lost.")
    print("  lambda* = the depth at which the floor starts to bind: one less")
    print("  the lowest (S - Em)/B over emissions in the unfloored run.\n")
    rows, out = [], {}
    dt = 0.02
    for name, Vf, prep, T in CASES:
        if q:
            T *= 0.5
        run = SeaLedger(Vf())
        e0 = prep(run)
        led = MinimalLedger(run)
        for _, force, absorb in MODES:
            lab = mode_label(force, absorb)
            base = led.run(e0, T, dt, sea_force=force, absorb=absorb)
            lam = base["deficit"]
            m0 = run.moments(base["e"])
            print(f"  {name}, {lab}: lambda* = {lam:.4f}, f = {base['f']:.4f},"
                  f" worst excess {base['excess']:.4f},"
                  f" |E-QLE| = {base['fid']:.3e}")
            print(f"  {'beta':>7} {'blocked':>9} {'eps':>9} {'dnorm':>9}"
                  f" {'d<p>':>9} {'dH':>9} {'dM-':>9} {'dT':>9}")
            betas = sorted({0.25, 0.5, 1.0, 1.5, round(0.98 * lam, 3),
                            round(1.02 * lam, 3), round(2.0 * lam, 3)})
            for beta in betas:
                o = led.run(e0, T, dt, sea_force=force, absorb=absorb,
                            beta=beta, rule="floor", ref=False)
                mm = run.moments(o["e"])
                eps = rel_l2(o["e"], base["e"])
                d = {k: mm[k] - m0[k] for k in mm}
                print(f"  {beta:7.3f} {o['blocked_frac']:9.2e} {eps:9.2e}"
                      f" {d['norm']:+9.1e} {d['p']:+9.1e} {d['H']:+9.1e}"
                      f" {d['Mneg']:+9.1e} {d['T']:+9.1e}")
                rows.append((name, lab, lam, beta, o["blocked_frac"], eps,
                             d["norm"], d["p"], d["H"], d["Mneg"], d["T"]))
                out[(name, lab, beta)] = eps
            out[(name, lab, "lam")] = lam
            print(flush=True)
    write_csv("sea_depletion_floor.csv",
              ["case", "mode", "lambda_star", "beta", "blocked_frac",
               "eps", "dnorm", "dp", "dH", "dMneg", "dT"], rows)
    print("  At beta >= lambda* the floor never binds and E is the unfloored")
    print("  E bit for bit.  Below it the lost emissions are lost")
    print("  depositions: the norm is untouched (each event deposits +1 and")
    print("  -1), while <p>, the energy, the negative mass and the")
    print("  transmission move.")
    return out


# ----------------------------------------------------------------------
# D.  Proposition V4: granular supply
# ----------------------------------------------------------------------
def part_d(q):
    banner("D  Proposition V4: granular supply at the physical density")
    print("  Q9(c): B dp 2 y_max = 1, so the event aperture holds one pair")
    print("  on average at the physical density.  An integer sea there is")
    print("  empty with probability exp(-n), n = nu S/B pairs.  Expected E")
    print("  under that rule against the sea depth beta, nu = 1.\n")
    rows, out = [], {}
    dt = 0.02
    nus = (1.0, 2.0, 4.0, 8.0, 16.0)
    betas = (1.0, 2.0, 4.0, 8.0, 16.0)
    for name, Vf, prep, T in CASES:
        if q:
            T *= 0.5
        run = SeaLedger(Vf())
        e0 = prep(run)
        led = MinimalLedger(run)
        for _, force, absorb in MODES[:2]:
            lab = mode_label(force, absorb)
            base = led.run(e0, T, dt, sea_force=force, absorb=absorb,
                           nus=nus, ref=False)
            est = base["est"]
            pb = est["blocked"] / max(est["demand"], 1e-30)
            print(f"  {name}, {lab}: emissive share 1 - f = "
                  f"{1 - base['f']:.4f}; share of emissions an integer sea"
                  " of depth 1 would block, from the unfloored run:")
            print("   " + "  ".join(f"nu={v:g}: {b:.3e}"
                                     for v, b in zip(nus, pb)))
            m0 = run.moments(base["e"])
            print(f"  {'beta':>7} {'blocked':>9} {'eps':>9} {'eps e^b':>9}"
                  f" {'d<p>':>9} {'dH':>9} {'dM-':>9} {'dT':>9}")
            for beta in betas:
                o = led.run(e0, T, dt, sea_force=force, absorb=absorb,
                            beta=beta, rule="poisson", nu=1.0, ref=False)
                mm = run.moments(o["e"])
                eps = rel_l2(o["e"], base["e"])
                print(f"  {beta:7.1f} {o['blocked_frac']:9.2e} {eps:9.2e}"
                      f" {eps * np.exp(beta):9.3f}"
                      f" {mm['p'] - m0['p']:+9.1e}"
                      f" {mm['H'] - m0['H']:+9.1e}"
                      f" {mm['Mneg'] - m0['Mneg']:+9.1e}"
                      f" {mm['T'] - m0['T']:+9.1e}")
                rows.append((name, lab, beta, o["blocked_frac"], eps,
                             mm["p"] - m0["p"], mm["H"] - m0["H"],
                             mm["Mneg"] - m0["Mneg"], mm["T"] - m0["T"]))
                out[(name, lab, beta)] = (eps, o["blocked_frac"])
            out[(name, lab, "est")] = pb
            print(flush=True)
    write_csv("sea_depletion_poisson.csv",
              ["case", "mode", "beta", "blocked_frac", "eps", "dp", "dH",
               "dMneg", "dT"], rows)
    return out


# ----------------------------------------------------------------------
# E.  lambda* over long runs and against the reach
# ----------------------------------------------------------------------
def part_e(q, t_long=24.0, reach_ladder=True):
    banner("E  lambda* over long runs (Q-SP2) and against the reach")
    rows, traces = [], {}
    T_long = 0.5 * t_long if q else t_long
    dt = 0.02
    print(f"  lambda*(t), the depth the floor needs so far, unfloored,"
          f" to T = {T_long}:")
    ts = [T_long * k / 6.0 for k in range(1, 7)]
    print(f"  {'case':>18} " + " ".join(f"{t:7.1f}" for t in ts)
          + f" {'excess':>8}")
    for name, Vf, prep, _ in CASES:
        run = SeaLedger(Vf())
        e0 = prep(run)
        led = MinimalLedger(run)
        for _, force, absorb in MODES:
            lab = mode_label(force, absorb)
            o = led.run(e0, T_long, dt, sea_force=force, absorb=absorb,
                        every=0.5, ref=False)
            tr = o["trace"]
            run_max = tr[:, 7]          # lambda*(t), cumulative
            traces[(name, lab)] = (tr[:, 0], tr[:, 4], run_max)
            vals = [run_max[np.argmin(np.abs(tr[:, 0] - t))] for t in ts]
            print(f"  {name[:7] + ' ' + lab:>18} "
                  + " ".join(f"{v:7.3f}" for v in vals)
                  + f" {o['excess']:8.3f}", flush=True)
            for t, v, x, rm in zip(tr[:, 0], tr[:, 4], tr[:, 5], run_max):
                rows.append((name, lab, t, v, x, rm))
    tag = "" if t_long == 24.0 else f"_T{t_long:g}"
    write_csv(f"sea_depletion_long{tag}.csv",
              ["case", "mode", "t", "deficit", "excess", "lambda_star"],
              rows)
    if not reach_ladder:
        return traces, []

    T = 4.0 if q else 8.0
    n_p, dpp = (64, 0.25) if q else (128, 0.125)
    yhs = (np.pi / 2, np.pi, 2 * np.pi) if q else (np.pi, 2 * np.pi,
                                                    4 * np.pi)
    print(f"\n  Against the reach: Eckart, minimal ledger, dp = {dpp},"
          f" T = {T}, window cut at y_h:")
    print(f"  {'y_h/a':>7} {'Gamma_tot max':>14} {'f (S′)':>8}"
          f" {'lam* (S)':>9} {'lam* (S′)':>10} {'lam* emi':>9}")
    reach = []
    for yh in yhs:
        run = SeaLedger(eckart(), n_p=n_p, dp=dpp, y_h=yh)
        e0 = summit_packet(run)
        led = MinimalLedger(run)
        res = [led.run(e0, T, dt, sea_force=f_, absorb=a_, ref=False)
               for _, f_, a_ in MODES]
        g = float(run.gamma_tot.max())
        print(f"  {yh:7.3f} {g:14.4f} {res[1]['f']:8.4f}"
              f" {res[0]['deficit']:9.4f} {res[1]['deficit']:10.4f}"
              f" {res[2]['deficit']:9.4f}", flush=True)
        reach.append((yh, g, res[1]["f"], res[0]["deficit"],
                      res[1]["deficit"], res[2]["deficit"]))
    write_csv("sea_depletion_reach.csv",
              ["y_h", "gamma_tot_max", "f_Sprime", "lambda_S",
               "lambda_Sprime", "lambda_emissive"], reach)
    return traces, reach


# ----------------------------------------------------------------------
# G.  Pair collisions: the amended (S′)  (addendum, Propositions V8-V10)
# ----------------------------------------------------------------------
class SeaMotion:
    """A W-null dynamics of the sea's pair density S(x, p), applied once per
    step after the events, on top of the (S′) streaming.  It moves whole
    pairs, so E never sees it (Proposition V8).

    kind 'none'      : (S′) as in step 23: rows keep their pairs.
         'volterra'  : linearised Volterra sea on the PAIR density,
                       dS/dt = K_res * S, the compensated residual acting on
                       S.  Antisymmetric, reversible.  Exact step in (x, s):
                       S_hat <- exp(h sym_e) S_hat.  (Proposition V9.)
         'diffusion' : symmetric control with the same channels and rates,
                       dS/dt = sum_q |K_q| (S_{n+q} + S_{n-q} - 2 S_n).
                       Exact step: S_hat <- exp(h (sym_a - Gamma_tot)) S_hat.
         'collide-k' : Definition (C), Focus/Defocus between aligned pairs,
                       mass action with detailed balance, on every channel q
                       at the kernel's rate: r_q B = scale |K_q(x)| / 2.
         'collide-u' : Definition (C) between neighbouring rows only (one
                       momentum step dp), at a rate independent of x and V:
                       r B = scale.
    Collisions per channel q, centre n:  J_n = r_q (S_n^2 - S_{n-q} S_{n+q})
    (defocus minus focus), dS_n/dt = -2 J_n + J_{n-q} + J_{n+q}.  Integrated
    by Heun's method with enough substeps for the linear symbol's largest
    rate, 16 r_q B summed over channels.
    """

    def __init__(self, run, kind="none", scale=1.0):
        self.kind, self.scale = kind, scale
        self.run = run
        if kind == "volterra":
            self.g_sym = run.sym_e
        elif kind == "diffusion":
            self.g_sym = run.sym_a - run.gamma_tot[:, None]
        elif kind == "collide-k":
            self.chan = [(q, scale * np.abs(run.k[:, q])[:, None] / (2.0 * B))
                         for q in range(1, run.n_p // 2)
                         if np.abs(run.k[:, q]).max() > 1e-14]
        elif kind == "collide-u":
            self.chan = [(1, np.full((run.r.size, 1), scale / B))]
        if kind.startswith("collide"):
            self.lam = float(sum(16.0 * r.max() * B for _, r in self.chan))

    def rhs(self, S):
        out = np.zeros_like(S)
        for q, r in self.chan:
            J = r * (S * S - np.roll(S, q, axis=1) * np.roll(S, -q, axis=1))
            out += -2.0 * J + np.roll(J, q, axis=1) + np.roll(J, -q, axis=1)
        return out

    def step(self, S, h):
        if self.kind == "none":
            return S
        if self.kind in ("volterra", "diffusion"):
            return np.real(np.fft.ifft(np.exp(h * self.g_sym)
                                       * np.fft.fft(S, axis=1), axis=1))
        n = max(1, int(np.ceil(h * self.lam / 0.5)))
        hs = h / n
        for _ in range(n):
            k1 = self.rhs(S)
            k2 = self.rhs(S + hs * k1)
            S = S + 0.5 * hs * (k1 + k2)
        return S


def sea_functionals(S, area):
    """Phi = sum [S - B - B ln(S/B)] (Volterra's invariant, Prop V9) and
    Hc = sum [S ln(S/B) - S + B] (the collisions' entropy, Prop V10), both
    in units of B per unit area; nan if any cell is empty or negative."""
    if S.min() <= 0.0:
        return float("nan"), float("nan")
    x = S / B
    return (float((x - 1.0 - np.log(x)).sum()) * area,
            float((x * np.log(x) - x + 1.0).sum()) * area)


def run_motion(run, e0, t_max, dt, motion, sea_force=False, absorb=True,
               every=0.5):
    """``MinimalLedger.run`` with the sea's own motion applied after the
    events of each step.  Returns the final E and a trace of
    (t, lambda*(t), worst cell 1 - min S/B, ||S - <S>_row||_2 / B, debt D,
    Phi, Hc, min S/B)."""
    run.sea_force = sea_force
    a = run.area
    E = e0.copy()
    S = np.full_like(e0, B)
    S0 = float(S.sum()) * a
    smin = [1.0]
    tr = []
    n_steps = int(round(t_max / dt))
    k_every = max(1, int(round(every / dt)))
    for step in range(n_steps + 1):
        if step % k_every == 0 or step == n_steps:
            dev = S - S.mean(axis=1, keepdims=True)
            Phi, Hc = sea_functionals(S, a)
            tr.append((step * dt, 1.0 - smin[0], 1.0 - float(S.min()) / B,
                       float(np.sqrt((dev ** 2).sum() * a)) / B,
                       S0 - float(S.sum()) * a, Phi, Hc, float(S.min()) / B))
        if step == n_steps:
            break
        for half in (0, 1):
            if half == 1:
                up, um = np.maximum(E, 0.0), np.maximum(-E, 0.0)
                up, um, S, *_ = run.channels_supply(up, um, S, dt, "none",
                                                    1.0, None, absorb,
                                                    smin=smin)
                c = np.minimum(up, um)
                S += c
                E = up - um
                S = motion.step(S, dt)
            E = run.stream(E, 0.5 * dt)
            S = run.stream_sea(S, 0.5 * dt)
    return E, np.array(tr)


def g_exact():
    """SymPy: the identities behind Propositions V9 and V10."""
    import sympy as sp
    print("  Exact (SymPy), periodic lattice of N = 7 rows, two channels:")
    N = 7
    S = sp.symbols("S0:%d" % N, positive=True)
    c1, c2, r1, r2, Bs, th = sp.symbols("c1 c2 r1 r2 B theta", positive=True)
    chans = ((1, c1), (2, c2))

    def volt(n):
        return S[n] / Bs * sum(c * (S[(n - q) % N] - S[(n + q) % N])
                               for q, c in chans)

    dsum = sp.simplify(sum(volt(n) for n in range(N)))
    dlog = sp.simplify(sum(volt(n) / S[n] for n in range(N)))
    print(f"    Volterra sea:  d/dt sum S = {dsum},  d/dt sum ln S = {dlog}"
          "  -> Phi = sum [S - B - B ln(S/B)] is conserved")

    def J(n, q, r):
        return r * (S[n % N] ** 2 - S[(n - q) % N] * S[(n + q) % N])

    def coll(n):
        return sum(-2 * J(n, q, r) + J(n - q, q, r) + J(n + q, q, r)
                   for q, r in ((1, r1), (2, r2)))

    msum = sp.expand(sum(coll(n) for n in range(N)))
    lhs = sum(sp.log(S[n]) * coll(n) for n in range(N))
    rhs = -sum(J(n, q, r) * (2 * sp.log(S[n]) - sp.log(S[(n - q) % N])
                             - sp.log(S[(n + q) % N]))
               for n in range(N) for q, r in ((1, r1), (2, r2)))
    hid = sp.simplify(sp.expand(lhs - rhs))
    print(f"    Collisions:    d/dt sum S = {msum};  dHc/dt + sum_n,q J (ln S_n^2"
          f" - ln S_n-q S_n+q) = {hid}")
    print("                   each term J (ln a - ln b) with J ∝ a - b is >= 0,"
          " so dHc/dt <= 0")
    n_, q_ = sp.symbols("n q", integer=True)
    print("    one Defocus event (two pairs at n -> n-q, n+q): change of"
          " sum n =", sp.simplify((n_ - q_) + (n_ + q_) - 2 * n_))
    # linearise at S = B + d: J ~ r B (2 d_n - d_{n-q} - d_{n+q}); symbol
    Jsym = r1 * Bs * (2 - 2 * sp.cos(th))
    sym = sp.simplify((-2 + 2 * sp.cos(th)) * Jsym)
    print(f"    linearised collision symbol per channel (theta = q x row"
          f" wavenumber): {sp.factor(sym)}  = -4 r B (1 - cos theta)^2")
    return dict(dsum=dsum, dlog=dlog, msum=msum, hid=hid)


def g_invariants():
    """Numerical: F conserved by the nonlinear Volterra sea, H decreasing
    under collisions, pairs and pair momentum conserved."""
    from scipy.integrate import solve_ivp
    rng = np.random.default_rng(7)
    n_p = 64
    S0 = B * (1.0 + 0.5 * rng.uniform(-1, 1, n_p))
    S0[10] = 0.05 * B
    c = rng.normal(size=6)

    def volt(t, S):
        acc = np.zeros_like(S)
        for q, cq in enumerate(c, start=1):
            acc += cq * (np.roll(S, q) - np.roll(S, -q))
        return S / B * acc
    sol = solve_ivp(volt, (0, 20), S0, rtol=1e-11, atol=1e-13,
                    method="DOP853", dense_output=True)
    Phi = lambda S: float(np.sum(S / B - 1 - np.log(S / B)))  # noqa: E731
    print("\n  Nonlinear Volterra sea, one row, six channels with random"
          " signed rates:")
    for t in (0, 5, 10, 20):
        S = sol.sol(t)
        print(f"    t = {t:4.0f}  sum S/B = {S.sum() / B:.10f}  Phi = {Phi(S):.10f}"
              f"  min S/B = {S.min() / B:.4f}")
    print("    -> Phi exact to the integrator's tolerance; the deep cell is"
          " shared out, Phi is not reduced.  An empty cell stays empty"
          " (dS_n/dt is proportional to S_n).")

    n_row = 1024                  # long enough that the ripple's tails
    #                               never reach the periodic row's ends,
    #                               where p jumps from +p_max to -p_max

    class _R:                     # one row, a minimal stand-in for SeaLedger
        n_p = n_row
        r = np.zeros(1)
        k = np.zeros((1, n_row))
    rr = _R()
    rr.k[0, 1:4] = (1.0, 0.6, 0.3)
    mo = SeaMotion(rr, "collide-k", 1.0)
    pgrid = np.arange(n_row) - n_row / 2
    S = np.full((1, n_row), B)
    S[0, 500:524] += B * 0.6 * np.sin(np.arange(24) / 3.0) * np.hanning(24)
    m0, p0 = S.sum(), (S * pgrid).sum()
    H0 = sea_functionals(S, 1.0)[1]
    print("\n  Collisions on one row of 1024 cells (channels 1-3), a ripple"
          " at its centre:")
    Hs = [H0]
    for k in range(1, 401):
        S = mo.step(S, 0.02)
        if k % 100 == 0:
            Hs.append(sea_functionals(S, 1.0)[1])
            print(f"    t = {k * 0.02:4.1f}  d(sum S)/sum S = "
                  f"{(S.sum() - m0) / m0:.1e}  d(sum p S) = "
                  f"{(S * pgrid).sum() - p0:.1e}  Hc = {Hs[-1]:.6f}"
                  f"  max|S/B - 1| = {np.abs(S / B - 1).max():.4f}")
    print(f"    -> Hc decreasing at every sample: {all(np.diff(Hs) < 0)}")


G_RUNS_MAIN = (("(S′) only", "none", 1.0), ("Volterra sea", "volterra", 1.0),
               ("pair diffusion", "diffusion", 1.0),
               ("collisions, kernel rate", "collide-k", 1.0),
               ("collisions, uniform rB = 1", "collide-u", 1.0))
G_SCAN = (("collide-k", (0.3, 0.1, 0.03, 0.01)),
          ("collide-u", (4.0, 0.25)))


def part_g(q, t_long=24.0, subs="IMSRX"):
    banner("G  Pair collisions: the amended (S′) (Propositions V8-V10)")
    T = 0.5 * t_long if q else t_long
    dt = 0.02
    tag = "" if t_long == 24.0 else f"_T{t_long:g}"
    if "I" in subs:
        g_exact()
        g_invariants()
    rows = []
    hdr = ["case", "realisation", "label", "kind", "scale", "t",
           "lambda_star", "worst", "l2dev", "debt", "Phi", "Hc", "smin",
           "max_dE"]

    def record(case, real, label, kind, scale, E, tr, E_ref):
        dE = (float(np.abs(E - E_ref).max()) if E_ref is not None
              else float("nan"))
        for row in tr:
            rows.append((case, real, label, kind, scale, *row, dE))
        return dE

    def flush(name):
        if rows:
            write_csv(f"sea_collisions_{name}{tag}.csv", hdr, rows)
        rows.clear()

    def show(label, tr, dE, ts):
        vals = [tr[np.argmin(np.abs(tr[:, 0] - t)), 1] for t in ts]
        print(f"  {label:>28} " + " ".join(f"{v:6.3f}" for v in vals)
              + f"  {tr[-1, 3]:7.3f} {tr[-1, 4]:6.3f}  {dE:7.1e}",
              flush=True)

    ts = [T * k / 6.0 for k in range(1, 7)]
    head = (f"  {'sea motion':>28} " + " ".join(f"{t:6.1f}" for t in ts)
            + f"  {'L2dev':>7} {'debt':>6}  {'max|dE|':>7}")

    def ladder(case, Vf, prep, absorb, runs, real):
        run = SeaLedger(Vf())
        e0 = prep(run)
        E_ref = None
        for label, kind, scale in runs:
            mo = SeaMotion(run, kind, scale)
            E, tr = run_motion(run, e0, T, dt, mo, False, absorb)
            if kind == "none":
                E_ref = E
            dE = record(case, real, label, kind, scale, E, tr, E_ref)
            show(label, tr, dE, ts)

    if "M" in subs:
        print(f"\n  lambda*(t), Pöschl–Teller well, minimal ledger, (S′),"
              f" absorptive first, to T = {T}; max|dE| against (S′) only:")
        print(head)
        ladder("PT well", well, well_packet, True, G_RUNS_MAIN, "abs")
        flush("main")
    if "S" in subs:
        print("\n  Rate scan, same run (rate in units of the kernel's,"
              " or rB for uniform):")
        print(head)
        runs = [("(S′) only", "none", 1.0)]
        for kind, scales in G_SCAN:
            runs += [(f"{kind} x {s:g}", kind, s) for s in scales]
        ladder("PT well", well, well_packet, True, runs, "abs")
        flush("scan")
    if "R" in subs:
        print("\n  Other cases, collisions at the kernel's rate:")
        print(head)
        rob = (("(S′) only", "none", 1.0),
               ("collisions, kernel rate", "collide-k", 1.0))
        for case, Vf, prep, absorb, real in (
                ("PT well", well, well_packet, False, "emi"),
                ("Eckart", eckart, summit_packet, True, "abs"),
                ("Eckart", eckart, summit_packet, False, "emi")):
            print(f"  -- {case}, {'absorptive first' if absorb else 'emissive'}")
            ladder(case, Vf, prep, absorb, rob, real)
        flush("cases")
    if "X" in subs:
        Tr = 4.0 if q else 8.0
        n_p, dpp = (64, 0.25) if q else (128, 0.125)
        yhs = (np.pi / 2, np.pi, 2 * np.pi) if q else (np.pi, 2 * np.pi,
                                                        4 * np.pi)
        print(f"\n  Against the reach: Eckart, (S′), dp = {dpp}, T = {Tr};"
              f" lambda* without and with collisions at the kernel's rate:")
        print(f"  {'y_h/a':>7} {'Gamma_tot max':>14} {'abs':>7} {'abs+C':>7}"
              f" {'emi':>7} {'emi+C':>7}")
        reach = []
        for yh in yhs:
            run = SeaLedger(eckart(), n_p=n_p, dp=dpp, y_h=yh)
            e0 = summit_packet(run)
            vals = []
            for absorb in (True, False):
                for kind in ("none", "collide-k"):
                    _, tr = run_motion(run, e0, Tr, dt,
                                       SeaMotion(run, kind, 1.0), False,
                                       absorb)
                    vals.append(tr[-1, 1])
            g = float(run.gamma_tot.max())
            print(f"  {yh:7.3f} {g:14.4f} " + " ".join(f"{v:7.3f}"
                                                       for v in vals),
                  flush=True)
            reach.append((yh, g, *vals))
        write_csv("sea_collisions_reach.csv",
                  ["y_h", "gamma_tot_max", "lambda_abs", "lambda_abs_C",
                   "lambda_emi", "lambda_emi_C"], reach)


def fig_collisions():
    """Drawn from the longest sea_collisions_{main,scan}*.csv present."""
    import glob
    import os
    best = []
    for part in ("main", "scan"):
        pick, t_best = None, -1.0
        for path in glob.glob(output_path(f"sea_collisions_{part}*.csv")):
            rows = read_csv(os.path.basename(path))
            t_max = max(float(r["t"]) for r in rows) if rows else -1.0
            if t_max > t_best:
                pick, t_best = rows, t_max
        best += pick or []
    if not best:
        return
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.2))
    pt = [r for r in best if r["case"] == "PT well" and r["realisation"]
          == "abs"]
    seen = set()
    cols = {"none": "k", "volterra": "C1", "diffusion": "C2",
            "collide-k": "C3", "collide-u": "C0"}
    for r in pt:
        key = r["label"]
        if key in seen:
            continue
        seen.add(key)
        rr = [x for x in pt if x["label"] == key]
        # a label can recur across sub-parts (the baseline): keep one run
        rr = rr[:len({x["t"] for x in rr})]
        t = [float(x["t"]) for x in rr]
        lam = [float(x["lambda_star"]) for x in rr]
        l2 = [float(x["l2dev"]) for x in rr]
        main = (r["label"] in [g[0] for g in G_RUNS_MAIN])
        style = dict(color=cols[r["kind"]], lw=1.6 if main else 0.9,
                     ls="-" if main else "--")
        lab = r["label"] if main else f"{r['kind']} × {float(r['scale']):g}"
        ax[0].plot(t, lam, label=lab if main else None, **style)
        ax[1].plot(t, l2, label=lab if main else None, **style)
        if not main:
            tag_ = (f"×{float(r['scale']):g}" if r["kind"] == "collide-k"
                    else f"rB={float(r['scale']):g}")
            ax[0].text(t[-1] + 1.0, lam[-1], tag_, fontsize=7,
                       color=style["color"], va="center")
            ax[1].text(t[-1] + 1.0, l2[-1], tag_, fontsize=7,
                       color=style["color"], va="center")
    ax[0].plot([], [], color=cols["collide-k"], ls="--", lw=0.9,
               label="rate scan, kernel rate × (end label)")
    ax[0].plot([], [], color=cols["collide-u"], ls="--", lw=0.9,
               label="rate scan, uniform (end label)")
    ax[0].axhline(1.0, color="k", lw=0.5, ls=":")
    ax[0].set_title("Pöschl–Teller well, (S′), absorptive first: λ*(t)")
    ax[1].set_title("sea deviation from its row means, ‖s − s̄‖₂ / B")
    for a_ in ax:
        a_.set_xlabel("t")
        a_.set_xlim(0, 106)
    ax[0].legend(fontsize=7)
    save_fig(fig, "sea_collisions.png")


# ----------------------------------------------------------------------
# R.  The transport's share of dM- against resolution (not in the default)
# ----------------------------------------------------------------------
def part_r(q):
    banner("R  The transport's share of dM- against resolution")
    print("  Eckart summit, minimal ledger, (S'), T = 8, reach 2 pi at every")
    print("  grid (window cut at y_h = 2 pi).\n")
    print(f"  {'n_r':>5} {'dr':>7} {'n_p':>5} {'dp':>7} {'D':>9} {'dM-':>9}"
          f" {'dN_tr/2':>9} {'|E-QLE|':>9}")
    grids = ((128, 64, 0.25), (256, 64, 0.25), (256, 128, 0.125),
             (512, 128, 0.125))
    for n_r, n_p, dp in (grids[:2] if q else grids):
        run = SeaLedger(eckart(), n_r=n_r, n_p=n_p, dp=dp, y_h=2 * np.pi)
        e0 = summit_packet(run)
        o = MinimalLedger(run).run(e0, 8.0, 0.02, sea_force=False,
                                   every=8.0)
        _, D, dm, ntr2 = o["trace"][-1][:4]
        print(f"  {n_r:5d} {run.dr:7.4f} {n_p:5d} {dp:7.4f} {D:9.5f}"
              f" {dm:9.5f} {ntr2:+9.5f} {o['fid']:9.2e}", flush=True)


# ----------------------------------------------------------------------
# F.  Demo defect: step 16's ledger books ringing as sea traffic
# ----------------------------------------------------------------------
def part_f(q):
    banner("F  Demo defect: transport ringing booked as sea traffic")
    print("  Step 16's ledger transports u+ and u- separately and")
    print("  spectrally.  Each is the positive part of E, kinked where E")
    print("  changes sign, so it rings, and its undershoots are negative.")
    print("  The allocation A = min(D, partners) then goes negative: it")
    print("  'absorbs' a negative partner, which adds a body pair and")
    print("  debits the sea with no event behind it.\n")
    T = np.pi if q else 2.0 * np.pi
    run = SeaLedger(harmonic())
    e0 = cat(run, 3.0, 1.0 / np.sqrt(2.0))
    o = run.run16(e0, T, 0.02)
    print(f"  cat in the harmonic trap, T = {T:.3f}, dt = 0.02 (no events:"
          f" max|K| = {np.abs(run.k[:, 1:]).max():.1e}):")
    print(f"    N(0) = {run.l1(e0):.4f}  N(T) = {o['N']:.4f}"
          f"   ringing pairs booked = {o['ring']:.4f}"
          f"   worst cell S/B = {o['min_s']:.4f}")
    o2 = MinimalLedger(run).run(e0, T, 0.02, ref=False)
    print(f"    minimal ledger: N(T) = {o2['N']:.4f}, worst cell S/B ="
          f" {float(o2['S'].min()) / B:.4f}, transport"
          f" dN_tr = {o2['n_tr']:+.4f}")
    print("\n  Step 16's Theorem S8 row (Eckart summit, dt = 0.01, T = 6,"
          " absorptive):")
    T6 = 3.0 if q else 6.0
    run = SeaLedger(eckart())
    e0 = summit_packet(run)
    print(f"  {'motion':>7} {'N(T)':>8} {'abs':>8} {'emi':>8} {'ringing':>9}"
          f" {'f':>7} {'min S/B':>8}")
    rows = []
    for force in (True, False):
        o = run.run16(e0, T6, 0.01, sea_force=force)
        f_ = o["n_abs"] / (o["n_abs"] + o["n_emi"])
        print(f"  {'(S)' if force else '(S′)':>7} {o['N']:8.4f}"
              f" {o['n_abs']:8.4f} {o['n_emi']:8.4f} {o['ring']:9.4f}"
              f" {f_:7.4f} {o['min_s']:8.4f}")
        rows.append(o)
        o2 = MinimalLedger(run).run(e0, T6, 0.01, sea_force=force,
                                    ref=False)
        print(f"  {'':>7} minimal ledger: N(T) = {o2['N']:.4f},"
              f" f = {o2['f']:.4f}, lambda* = {o2['deficit']:.4f}")
    print("\n  The 'abs' column counts the allocation with its sign; the")
    print("  ringing column is the negative part inside it.  E is not")
    print("  affected (A + Em = D whatever A is), but N, S and f are.")
    return rows


# ----------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------
def fig_free(a_out):
    names = ("cat, harmonic", "cat, PT well")
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    run, e0, _ = a_out[names[0]]
    sl = np.argsort(run.p)
    ext = (run.p[sl][0], run.p[sl][-1], run.r[0], run.r[-1])
    vm = float(np.abs(e0).max())
    ax[0].imshow(e0[:, sl], origin="lower", aspect="auto", extent=ext,
                 cmap="RdBu_r", vmin=-vm, vmax=vm)
    ax[0].set_title("even cat at t = 0")
    for i, nm in enumerate(names):
        run, e0, o = a_out[nm]
        ax[i + 1].imshow(o["e"][:, sl], origin="lower", aspect="auto",
                         extent=ext, cmap="RdBu_r", vmin=-vm, vmax=vm)
        ax[i + 1].set_title(f"{nm}, t = 2π: events {o['n_ev']:.2g}")
    for a in ax:
        a.set_xlabel("p")
        a.set_ylabel("x")
        a.set_ylim(-7, 7)
        a.set_xlim(-5, 5)
    fig.suptitle("Proposition V1: the trap spends no sea pair on the cat's "
                 "negativity; the anharmonic well spends many", y=1.03)
    save_fig(fig, "sea_depletion_free_negativity.png")


def read_csv(name):
    import os
    path = output_path(name)
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return list(csv.DictReader(fh))


MODE_STYLE = {"(S) abs": "o-", "(S′) abs": "s--", "(S′) emi": "^:"}
CASE_COL = {"Eckart": "C0", "PT well": "C3"}


def fig_floor():
    """Built from the CSVs of Parts C and D, so they can run separately."""
    fc, fd = read_csv("sea_depletion_floor.csv"), read_csv(
        "sea_depletion_poisson.csv")
    if fc is None:
        return
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.4))
    for name, col in CASE_COL.items():
        for lab, mk in MODE_STYLE.items():
            rr = [r for r in fc if r["case"] == name and r["mode"] == lab]
            if not rr:
                continue
            bs = [float(r["beta"]) for r in rr]
            ep = [max(float(r["eps"]), 1e-17) for r in rr]
            ax[0].semilogy(bs, ep, mk, color=col, ms=4,
                           label=f"{name} {lab}")
            ax[0].axvline(float(rr[0]["lambda_star"]), color=col, lw=0.6,
                          ls=mk[1:] or "-")
    ax[0].set_xlabel("sea depth β (units of B)")
    ax[0].set_ylabel("ε = |E − E_unfloored| / |E_unfloored|")
    ax[0].set_title("V3: floor, mean field (ε = 0 drawn at 1e-17;"
                    " verticals: λ*)")
    ax[0].legend(fontsize=7)
    if fd is not None:
        for name, col in CASE_COL.items():
            for lab, mk in MODE_STYLE.items():
                rr = [r for r in fd if r["case"] == name
                      and r["mode"] == lab]
                if not rr:
                    continue
                ax[1].semilogy([float(r["beta"]) for r in rr],
                               [float(r["eps"]) for r in rr], mk, color=col,
                               ms=4, label=f"{name} {lab}")
        bb = np.linspace(1, 16, 50)
        ax[1].semilogy(bb, np.exp(-bb), "k:", label="exp(−β)")
        ax[1].set_xlabel("sea depth β (pairs per event aperture)")
        ax[1].set_title("V4: integer sea, expected E")
        ax[1].legend(fontsize=7)
    save_fig(fig, "sea_depletion_floor.png")


def fig_long():
    """Drawn from the longest long-run CSV present."""
    import glob
    import os
    best, t_best = None, -1.0
    for path in glob.glob(output_path("sea_depletion_long*.csv")):
        rows = read_csv(os.path.basename(path))
        t_max = max(float(r["t"]) for r in rows) if rows else -1.0
        if t_max > t_best:
            best, t_best = rows, t_max
    fl = best
    if fl is None:
        return
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4))
    for i, name in enumerate(CASE_COL):
        for lab, col in zip(MODE_STYLE, ("C0", "C3", "C2")):
            rr = [r for r in fl if r["case"] == name and r["mode"] == lab]
            if not rr:
                continue
            t = np.array([float(r["t"]) for r in rr])
            d = np.array([float(r["deficit"]) for r in rr])
            lam = np.array([float(r["lambda_star"]) for r in rr])
            ax[i].plot(t, d, color=col, lw=0.7, alpha=0.6)
            ax[i].plot(t, lam, color=col, lw=1.6, label=lab)
        ax[i].axhline(1.0, color="k", lw=0.5, ls=":")
        ax[i].set_title(f"{name}: worst cell 1 − S/B after each step"
                        " (thin) and λ*(t) (thick)")
        ax[i].set_xlabel("t")
        ax[i].legend(fontsize=8)
    save_fig(fig, "sea_depletion_long.png")


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--parts", default="ABCDEFG",
                    help="parts to run; '' redraws the CSV figures only")
    ap.add_argument("--quick", action="store_true",
                    help="halve run lengths (testing only)")
    ap.add_argument("--t-long", type=float, default=24.0,
                    help="Part E: length of the long runs (the note also"
                         " quotes 96)")
    ap.add_argument("--no-reach", action="store_true",
                    help="Part E: skip the reach ladder")
    ap.add_argument("--g-subs", default="IMSRX",
                    help="Part G sub-parts: I exact and invariants, M main"
                         " comparison, S rate scan, R other cases, X reach"
                         " (its long runs use --t-long)")
    args = ap.parse_args()
    t0 = time.time()
    if "A" in args.parts:
        fig_free(part_a(args.quick))
    if "B" in args.parts:
        part_b(args.quick)
    if "C" in args.parts:
        part_c(args.quick)
    if "D" in args.parts:
        part_d(args.quick)
    if "E" in args.parts:
        part_e(args.quick, args.t_long, not args.no_reach)
    if "F" in args.parts:
        part_f(args.quick)
    if "G" in args.parts:
        part_g(args.quick, args.t_long, args.g_subs)
    if "R" in args.parts:
        part_r(args.quick)
    fig_floor()                 # from whatever CSVs exist so far
    fig_long()
    fig_collisions()
    print(f"\nDone in {time.time() - t0:.0f} s.")


if __name__ == "__main__":
    main()
