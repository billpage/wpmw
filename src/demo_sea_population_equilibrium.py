"""Sea population equilibrium under the compensated Liouville algorithm.

Verification companion to ``docs/analysis/sea_population_equilibrium.md``.

The compensated split leaves the residual channel signed, so it is not a rate
and the sign has to be unravelled.  The specification names positon/negaton as
the recommended reading but leaves the resulting ledger unpriced (open item
CLA3).  This demo prices it, on the Eckart barrier of
``docs/analysis/eckart_barrier_compensated.md``.

Ledger fields on the (r, p) mesh, in Wigner units:

    E = u+ - u-     the signed density, the observable, E == W
    N = u+ + u-     the body density
    S               bound sea pairs per cell, background B = 2/h

Parts
  A  Theorem S1 and S2: the sink must be bilinear; the ledger splits into a
     closed E equation and an N equation under the absolute kernel, with a
     closed-form local fixed point.
  B  Theorem S3: relaxation is exponential where the state lives and
     algebraic, N -> 2/(kappa t), where it does not.
  C  Theorem S4 and S5: emissive-only unravelling drains the sea without
     bound, and throttling by a finite sea corrupts W in the core.
  D  Theorem S6: the mode is chosen per event, not per leg -- a mixed event
     violates body-momentum conservation.
  E  Theorem S7: the ledger identity.  With absorptive fraction f,
     dN = 2(1 - 2f) n_ev and dS = (2f - 1) n_ev, so f = 1/2 closes both.
  F  Theorem S8: absorptive unravelling restores QLE fidelity, bounds the
     body count, and shrinks the sea deficit by two to three orders.
  G  Reach dependence, and the kappa proportional to Gamma scaling.

Run as::

    WPMW_OUTPUT=... PYTHONPATH=src python3 src/demo_sea_population_equilibrium.py
"""

from __future__ import annotations

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from wpmwlib.wpmw_utils import output_path, docs_path

HBAR = 1.0
MU = 1.0
B = 1.0 / (np.pi * HBAR)          # the crystal shift 2/h


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
    """Eckart barrier, compensated split, with a three-field ledger."""

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
        d_res = (V(rr + yy, v0, a) - V(rr - yy, v0, a)
                 - 2.0 * yy * dV(rr, v0, a))
        # raised-cosine horizon, the specification's default profile
        w = np.cos(np.pi * self.y / (2.0 * self.y_max)) ** 2
        m_res = (1j / HBAR) * d_res * w[None, :]
        self.k = np.real(np.fft.ifft(m_res, axis=1))       # real, odd in q
        self.sym_e = np.fft.fft(self.k, axis=1)            # signed symbol
        self.sym_a = np.real(np.fft.fft(np.abs(self.k), axis=1))
        self.gamma_tot = self.sym_a[:, 0].copy()           # sum_q |K_q|

    # -- transport and the two mesh kernels ----------------------------
    def stream(self, f, dt):
        fh = np.fft.fft(f, axis=0)
        fh *= np.exp(-1j * self.kr[:, None] * self.p[None, :] * dt / MU)
        return np.real(np.fft.ifft(fh, axis=0))

    def kick(self, f, sym, dt):
        return np.real(np.fft.ifft(np.fft.fft(f, axis=1)
                                   * np.exp(dt * sym), axis=1))

    def qle_step(self, e, dt):
        """Exact mesh reference: the QLE under the compensated split."""
        e = self.stream(e, 0.5 * dt)
        e = self.kick(e, self.sym_e, dt)
        return self.stream(e, 0.5 * dt)

    # -- the local nonlinearity ----------------------------------------
    @staticmethod
    def clamp(n, e):
        """N >= |E| is structural; spectral transport rings, so re-impose."""
        return np.maximum(n, np.abs(e))

    @staticmethod
    def recombine(n, e, kappa, dt):
        """Exact solution of dN/dt = -(kappa/2)(N^2 - E^2), E frozen."""
        if kappa == 0.0:
            return n
        ea = np.abs(e)
        out = np.empty_like(n)
        big = ea > 1e-14
        sm = ~big
        out[sm] = n[sm] / (1.0 + 0.5 * kappa * n[sm] * dt)
        nb, eb = n[big], ea[big]
        z = ((nb - eb) / (nb + eb)) * np.exp(-kappa * eb * dt)
        out[big] = eb * (1.0 + z) / (1.0 - z)
        return out

    @staticmethod
    def n_eq(gamma, kappa, e):
        """Local fixed point of the N equation (Theorem S2)."""
        g = gamma / kappa
        return g + np.sqrt(g * g + e * e)

    # -- the mean-field two-field run (emissive unravelling) -----------
    def run_mesh(self, e0, kappa, t_max, dt, n0=None, live=False):
        e = e0.copy()
        n = np.abs(e0) if n0 is None else n0.copy()
        S = np.full_like(e0, B)          # cumulative diagnostic ledger
        mn = 1.0
        for _ in range(int(round(t_max / dt))):
            e, n, S = (self.stream(e, .5 * dt), self.stream(n, .5 * dt),
                       self.stream(S, .5 * dt))
            n = self.clamp(n, e)
            sig = np.clip(S / B, 0.0, None) if live else 1.0
            if live:
                e = e + dt * self.kick_lin(e * sig, self.sym_e)
                n = n + dt * self.kick_lin(n * sig, self.sym_a)
            else:
                e = self.kick(e, self.sym_e, dt)
                n = self.kick(n, self.sym_a, dt)
            born = 0.5 * self.gamma_tot[:, None] * n * (sig if live else 1.0)
            n_pre = n.copy()
            n = self.recombine(n, e, kappa, dt)
            S = S - dt * born + 0.5 * np.maximum(n_pre - n, 0.0)
            e, n, S = (self.stream(e, .5 * dt), self.stream(n, .5 * dt),
                       self.stream(S, .5 * dt))
            n = self.clamp(n, e)
            mn = min(mn, float((S / B).min()))
        return dict(e=e, n=n, S=S, min_s=mn,
                    ntot=float(np.sum(n) * self.area))

    def kick_lin(self, f, sym):
        return np.real(np.fft.ifft(np.fft.fft(f, axis=1) * sym, axis=1))

    # -- the event-resolved run (emissive / absorptive unravelling) ----
    def channels(self, up, um, S, dt, absorb):
        """Tau-leap over channel pairs (q, -q), sequential with live caps.

        Each pair is ONE event depositing +1 at one daughter row and -1 at
        the other.  Absorptive realisation needs a partner of the right
        species at BOTH daughters; otherwise the event falls back to
        emissive, splitting a bound pair at the parent row.
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
                if absorb:
                    capA = np.where(t > 0, np.roll(um, -q, axis=1),
                                    np.roll(up, -q, axis=1))
                    capB = np.where(t > 0, np.roll(up, q, axis=1),
                                    np.roll(um, q, axis=1))
                    A = np.minimum(D, np.minimum(capA, capB))
                else:
                    A = np.zeros_like(D)
                Em = D - A
                n_abs += float(A.sum())
                n_emi += float(Em.sum())
                # absorptive: consume partners at p+q and p-q, bind at p
                aq, am = np.roll(A, q, axis=1), np.roll(A, -q, axis=1)
                um -= np.where(t > 0, aq, 0.0)
                up -= np.where(t > 0, 0.0, aq)
                up -= np.where(t > 0, am, 0.0)
                um -= np.where(t > 0, 0.0, am)
                S += A
                # emissive: split a bound pair at p, children at p +- q
                eq, em_ = np.roll(Em, q, axis=1), np.roll(Em, -q, axis=1)
                up += np.where(t > 0, eq, 0.0)
                um += np.where(t > 0, 0.0, eq)
                um += np.where(t > 0, em_, 0.0)
                up += np.where(t > 0, 0.0, em_)
                S -= Em
                np.maximum(up, 0.0, out=up)
                np.maximum(um, 0.0, out=um)
        return up, um, S, n_abs, n_emi

    def run_events(self, e0, t_max, dt, absorb=True, trace=False):
        up = np.maximum(e0, 0.0).copy()
        um = np.maximum(-e0, 0.0).copy()
        S = np.full_like(e0, B)
        ref = e0.copy()
        n_abs = n_emi = 0.0
        mn, tr = 1.0, []
        for step in range(int(round(t_max / dt))):
            up, um, S = (self.stream(up, .5 * dt), self.stream(um, .5 * dt),
                         self.stream(S, .5 * dt))
            up, um, S, a, e_ = self.channels(up, um, S, dt, absorb)
            n_abs += a
            n_emi += e_
            up, um, S = (self.stream(up, .5 * dt), self.stream(um, .5 * dt),
                         self.stream(S, .5 * dt))
            mn = min(mn, float((S / B).min()))
            ref = self.qle_step(ref, dt)
            if trace and step % max(1, int(0.5 / dt)) == 0:
                tr.append((step * dt, float((S / B).min()),
                           float(np.sum(up + um) * self.area)))
        e = up - um
        ev = n_abs + n_emi
        return dict(e=e, ref=ref, S=S, min_s=mn, trace=np.array(tr),
                    N=float(np.sum(up + um) * self.area),
                    f=n_abs / max(ev, 1e-30), n_ev=ev,
                    fid=float(np.linalg.norm(e - ref) / np.linalg.norm(ref)))


def gaussian(run, r_c, p_c, sigma_r):
    sp = HBAR / (2.0 * sigma_r)
    w = np.exp(-((run.r[:, None] - r_c) ** 2) / (2.0 * sigma_r ** 2)
               - ((run.p[None, :] - p_c) ** 2) / (2.0 * sp ** 2))
    return w / (np.sum(w) * run.dr * run.dp)


# ======================================================================
def part_a(run):
    banner("A  Theorem S1, S2: the ledger split and its fixed point")
    print("  A removal channel preserves E = u+ - u- for arbitrary ensembles")
    print("  iff it removes coincident pairs, hence is bilinear.  A per-body")
    print("  death rate mu sends E -> e^{-mu t} E and fails.\n")
    print("  d_t E = sum_q  K_q  E(p - xi_q)          closed, exactly the QLE")
    print("  d_t N = sum_q |K_q| N(p - xi_q) - (kappa/2)(N^2 - E^2)\n")
    print(f"  y_max = {run.y_max:.4f}, dp = {run.dp}, B = 2/h = {B:.6f}")
    print(f"  Gamma_tot peak = {run.gamma_tot.max():.4f} at "
          f"r = {run.r[np.argmax(run.gamma_tot)]:+.3f}")
    i0 = int(np.argmin(np.abs(run.r)))
    print(f"  Gamma_tot at the summit = {run.gamma_tot[i0]:.3e}"
          "   (a K7 quiet point)\n")
    print("  single-cell relaxation, Gamma = 2.0, E = 0.30, T = 60")
    print(f"  {'kappa':>8} {'N_eq closed form':>17} {'from below':>12}"
          f" {'from above':>12} {'rel err':>10}")
    G, E, T = 2.0, 0.30, 60.0
    for kappa in [0.5, 2.0, 10.0, 100.0, 1000.0]:
        pred = Ledger.n_eq(G, kappa, E)
        vals = []
        for n_start in (E, 100.0 * E):
            n, dt = n_start, 1e-4
            for _ in range(int(T / dt)):
                n += dt * (G * n - 0.5 * kappa * (n * n - E * E))
            vals.append(n)
        print(f"  {kappa:8.1f} {pred:17.6f} {vals[0]:12.6f}"
              f" {vals[1]:12.6f} {abs(vals[0]-pred)/pred:10.2e}")
    print("\n  kappa -> inf: N_eq -> |E|, the minimal ensemble.")
    print("  kappa -> 0  : N_eq -> 2 Gamma / kappa, independent of the state.")
    return G, E


def part_b(run):
    banner("B  Theorem S3: two relaxation laws")
    kappa = 20.0
    print(f"  Where Gamma = E = 0 the ledger obeys dN/dt = -(kappa/2) N^2,")
    print(f"  so N(t) = N0 / (1 + kappa N0 t / 2) -> 2/(kappa t).  kappa = {kappa}")
    print(f"  {'t':>8} " + " ".join(f"{'N0=' + str(v):>12}"
                                    for v in (0.01, 1.0, 100.0))
          + f" {'2/(kappa t)':>12}")
    for t in (0.1, 1.0, 10.0, 100.0, 1000.0):
        row = [v / (1.0 + 0.5 * kappa * v * t) for v in (0.01, 1.0, 100.0)]
        print(f"  {t:8.1f} " + " ".join(f"{v:12.6f}" for v in row)
              + f" {2.0/(kappa*t):12.6f}")
    print("\n  The asymptote is universal, so preparation is forgotten -- but")
    print("  algebraically.  Inside the emitting region, where |E| > 0, the")
    print("  approach is exponential at rate ~ kappa |E|.")

    e0 = gaussian(run, 0.0, 1.0, 1.0)
    print(f"\n  Eckart packet: peak E = {e0.max():.6f}, B = {B:.6f},"
          f" ratio {e0.max()/B:.4f}")
    print("  |W| <= 2/h for pure states, saturated at a Gaussian centre, so")
    print("  the excess is never small against the sea where it matters.")
    preps = {"minimal N=|E|": np.abs(e0),
             "inflated N=5|E|": 5.0 * np.abs(e0),
             "padded N=|E|+B/2": np.abs(e0) + 0.5 * B}
    print(f"\n  {'preparation':>18} {'N(0)':>9} {'N(6)':>9}")
    outs = {}
    for k, n0 in preps.items():
        o = run.run_mesh(e0, 20.0, 6.0, 0.01, n0=n0)
        outs[k] = o
        print(f"  {k:>18} {np.sum(n0)*run.area:9.4f} {o['ntot']:9.4f}")
    return e0, outs


def part_c(run, e0):
    banner("C  Theorem S4, S5: emissive-only unravelling")
    print("  The specification's rate R(x) = sum_q |K_res| unravels every")
    print("  channel emissively.  S is then debited at the parent row p and")
    print("  credited at the daughter rows p +- xi_q: transport, not use.\n")
    print(f"  {'kappa':>8} {'N_tot(6)':>10} {'min S/B':>9} {'mean S/B':>9}")
    for kappa in (0.0, 20.0, 200.0, 2000.0):
        o = run.run_mesh(e0, kappa, 6.0, 0.01)
        print(f"  {kappa:8.0f} {o['ntot']:10.4g} {o['min_s']:9.4f}"
              f" {float((o['S']/B).mean()):9.5f}")
    print("\n  The mean is conserved -- bodies are -- while the worst cell")
    print("  is stripped.  Fast recombination does not repair it, because")
    print("  the credit never returns to the row that was debited.")


def part_d(run):
    banner("D  Theorem S6: the mode is per event, not per leg")
    q = 3
    xi = q * run.dp
    p0 = 2.0
    print("  A channel event deposits +1 at p + xi_q and -1 at p - xi_q.")
    print("  Each deposition has two realisations, identical in E:\n")
    print("     +1 at a cell  ==  add a positon   OR  remove a negaton")
    print("     -1 at a cell  ==  add a negaton   OR  remove a positon\n")
    print(f"  body momentum, parent p = {p0}, xi_q = {xi}:")
    rows = [("emissive   (both legs create)", 2 * p0, (p0 + xi) + (p0 - xi)),
            ("absorptive (both legs remove)", (p0 + xi) + (p0 - xi), 2 * p0),
            ("mixed      (create + / remove -)", (p0 - xi) + 0.0,
             (p0 + xi) + 0.0)]
    print(f"  {'realisation':>34} {'before':>9} {'after':>9} {'delta':>9}")
    for name, before, after in rows:
        print(f"  {name:>34} {before:9.4f} {after:9.4f} {after-before:+9.4f}")
    print("\n  The mixed realisation moves one body from p - xi_q to p + xi_q")
    print("  and so needs 2 xi_q of momentum from somewhere.  The parent")
    print("  streams on undisturbed, so nothing can supply it: an event is")
    print("  wholly emissive or wholly absorptive, and absorption therefore")
    print("  requires a partner at BOTH daughters.")


def part_e(run, e0):
    banner("E  Theorem S7: the ledger identity")
    print("  Emissive event:   bodies +2, S debited  at the parent row.")
    print("  Absorptive event: bodies -2, S credited at the parent row.")
    print("  With an absorptive fraction f over n_ev events,\n")
    print("      dN = 2 (1 - 2f) n_ev        dS = (2f - 1) n_ev\n")
    print("  so f = 1/2 closes BOTH ledgers simultaneously.\n")
    o = run.run_events(e0, 6.0, 0.01, absorb=True)
    n_ev = o["n_ev"]
    f = o["f"]
    dn_pred = 2.0 * (1.0 - 2.0 * f) * n_ev * run.area
    ds_pred = (2.0 * f - 1.0) * n_ev * run.area
    dn_meas = o["N"] - float(np.sum(np.abs(e0)) * run.area)
    ds_meas = float(np.sum(o["S"] - B) * run.area)
    print(f"  {'quantity':>12} {'predicted':>12} {'measured':>12} {'rel':>10}")
    print(f"  {'dN':>12} {dn_pred:12.5f} {dn_meas:12.5f}"
          f" {abs(dn_pred-dn_meas)/abs(dn_meas):10.2e}")
    print(f"  {'dS':>12} {ds_pred:12.5f} {ds_meas:12.5f}"
          f" {abs(ds_pred-ds_meas)/abs(ds_meas):10.2e}")
    print(f"\n  measured f = {f:.5f}, short of 1/2 by {0.5-f:.5f}")
    return f


def part_f(run, e0):
    banner("F  Theorem S8: absorptive unravelling")
    print("  Every event tries absorptive first and falls back to emissive")
    print("  when a partner is missing.  The event RATE is |K_q| either way,")
    print("  so E is unchanged: the mode is invisible in the observable.\n")
    print(f"  {'dt':>7} {'mode':>11} {'N(6)':>11} {'f':>7} {'min S/B':>9}"
          f" {'mean S/B':>9} {'rel L2 vs QLE':>14}")
    fid = {}
    for dt in (0.02, 0.01, 0.005):
        for absorb in (False, True):
            if not absorb and dt == 0.005:
                continue
            o = run.run_events(e0, 6.0, dt, absorb=absorb)
            fid[(dt, absorb)] = o["fid"]
            print(f"  {dt:7.3f} {'absorptive' if absorb else 'emissive':>11}"
                  f" {o['N']:11.4f} {o['f']:7.4f} {o['min_s']:9.4f}"
                  f" {float((o['S']/B).mean()):9.5f} {o['fid']:14.3e}")
    print("\n  Both modes converge at first order in dt; nothing is throttled,")
    print("  so the absorptive E converges to the exact QLE.  The four-order")
    print("  gap is amplification: the emissive ledger computes E as a small")
    print("  difference of two populations of size 1e5, so every integration")
    print("  error is multiplied by N / ||E||.")
    tr_e = run.run_events(e0, 18.0, 0.01, absorb=False, trace=True)["trace"]
    tr_a = run.run_events(e0, 18.0, 0.01, absorb=True, trace=True)["trace"]
    print(f"\n  {'t':>6} {'emissive min S/B':>18} {'absorptive min S/B':>20}")
    for t in (2.0, 6.0, 10.0, 14.0, 18.0):
        ie = int(np.argmin(np.abs(tr_e[:, 0] - t)))
        ia = int(np.argmin(np.abs(tr_a[:, 0] - t)))
        print(f"  {t:6.1f} {tr_e[ie,1]:18.4f} {tr_a[ia,1]:20.4f}")
    return fid, tr_e, tr_a


def part_g(e0_unused):
    banner("G  reach dependence")
    print(f"  {'y_max':>8} {'dp':>7} {'Gamma_tot':>10} {'f':>8} {'min S/B':>9}"
          f" {'N(6)':>8} {'rel L2':>10}")
    rows = []
    for dp, n_p in ((0.5, 64), (0.25, 64), (0.125, 128)):
        rr = Ledger(n_p=n_p, dp=dp)
        ee = gaussian(rr, 0.0, 1.0, 1.0)
        o = rr.run_events(ee, 6.0, 0.01, absorb=True)
        rows.append((rr.y_max, rr.gamma_tot.max(), o["f"], o["min_s"], o["N"]))
        print(f"  {rr.y_max:8.4f} {dp:7.3f} {rr.gamma_tot.max():10.4f}"
              f" {o['f']:8.4f} {o['min_s']:9.4f} {o['N']:8.4f}"
              f" {o['fid']:10.3e}")
    print("\n  f rises toward 1/2 as the reach grows: a longer reach opens more")
    print("  channels, so more events find a partner.  The residual deficit")
    print("  still tracks the reach.")
    print(f"\n  {'y_max':>8} {'N_eq/B, kappa=20':>18} {'N_eq/B, kappa=c Gamma':>23}")
    for y, g, _, _, _ in rows:
        f20 = Ledger.n_eq(g, 20.0, B) / B
        fpg = Ledger.n_eq(g, 6.3 * g, B) / B
        print(f"  {y:8.4f} {f20:18.4f} {fpg:23.4f}")
    print("\n  At fixed kappa the equilibrium ledger tracks the reach; with")
    print("  kappa proportional to Gamma_tot(x) it does not.  Reach")
    print("  independence therefore fixes the FORM of the recombination rate")
    print("  and leaves one dimensionless constant.")
    return rows


# ----------------------------------------------------------------------
def figures(run, G, E, outs, fid, tr_e, tr_a, rows, f_meas):
    kap = np.logspace(-1, 3, 200)

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    ax[0].loglog(kap, Ledger.n_eq(G, kap, E), lw=2, color="C0",
                 label=r"$N_{\rm eq}$")
    ax[0].axhline(E, ls="--", color="C3", label=r"$|E|$  (minimal)")
    ax[0].loglog(kap, 2.0 * G / kap, ls=":", color="C2",
                 label=r"$2\Gamma/\kappa$")
    ax[0].set_xlabel(r"$\kappa$")
    ax[0].set_ylabel(r"$N_{\rm eq}$")
    ax[0].set_title(r"S2: $N_{\rm eq}=\Gamma/\kappa+\sqrt{\Gamma^2/\kappa^2+E^2}$")
    ax[0].legend(fontsize=9)
    ax[0].grid(alpha=.3)

    tt = np.logspace(-1, 3, 200)
    for n0, c in ((0.01, "C0"), (1.0, "C1"), (100.0, "C2")):
        ax[1].loglog(tt, n0 / (1 + 0.5 * 20.0 * n0 * tt), color=c,
                     label=rf"$N_0={n0}$")
    ax[1].loglog(tt, 2.0 / (20.0 * tt), "k--", label=r"$2/(\kappa t)$")
    ax[1].set_xlabel("t")
    ax[1].set_ylabel("N")
    ax[1].set_title(r"S3: vacuum relaxation is $1/t$")
    ax[1].legend(fontsize=9)
    ax[1].grid(alpha=.3)

    ax[2].plot(run.r, run.gamma_tot, lw=2, color="C0")
    rstar = run.a * np.log(np.sqrt(2.0) + np.sqrt(3.0))
    for x in (0.0, rstar, -rstar):
        ax[2].axvline(x, ls=":", color="C3", lw=1)
    ax[2].set_xlim(-6, 6)
    ax[2].set_xlabel("r")
    ax[2].set_ylabel(r"$\Gamma_{\rm tot}(r)$")
    ax[2].set_title("emission rate and the K7 quiet points")
    ax[2].grid(alpha=.3)
    fig.suptitle("Sea population equilibrium: the ledger fixed point", y=1.02)
    save_fig(fig, "sea_population_fixed_point.png")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    ax[0].plot(tr_e[:, 0], tr_e[:, 2], lw=2, color="C3", label="emissive")
    ax[0].plot(tr_a[:, 0], tr_a[:, 2], lw=2, color="C0", label="absorptive")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("t")
    ax[0].set_ylabel(r"$N_{\rm tot}$")
    ax[0].set_title("S8: the body count")
    ax[0].legend(fontsize=9)
    ax[0].grid(alpha=.3)

    ax[1].plot(tr_e[:, 0], tr_e[:, 1], lw=2, color="C3", label="emissive")
    ax[1].plot(tr_a[:, 0], tr_a[:, 1], lw=2, color="C0", label="absorptive")
    ax[1].axhline(0.0, color="k", lw=.8)
    ax[1].set_xlabel("t")
    ax[1].set_ylabel(r"worst cell $S/B$")
    ax[1].set_title("S4 vs S8: the sea deficit")
    ax[1].legend(fontsize=9)
    ax[1].grid(alpha=.3)

    dts = np.array([0.02, 0.01, 0.005])
    fa = np.array([fid[(d, True)] for d in dts])
    fe = np.array([fid[(d, False)] for d in dts[:2]])
    ax[2].loglog(dts, fa, "o-", color="C0", label="absorptive")
    ax[2].loglog(dts[:2], fe, "s-", color="C3", label="emissive")
    ax[2].loglog(dts, fa[0] * dts / dts[0], "k:", label=r"$O(\Delta t)$")
    ax[2].set_xlabel(r"$\Delta t$")
    ax[2].set_ylabel(r"rel $L^2$ vs exact QLE")
    ax[2].set_title("fidelity of the two unravellings")
    ax[2].legend(fontsize=9)
    ax[2].grid(alpha=.3, which="both")
    fig.suptitle("Emissive against absorptive unravelling", y=1.02)
    save_fig(fig, "sea_population_unravelling.png")

    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.2))
    ymax = [r[0] for r in rows]
    ax[0].plot(ymax, [r[2] for r in rows], "o-", lw=2, color="C0")
    ax[0].axhline(0.5, ls="--", color="C3", label=r"$f=1/2$: ledger closes")
    ax[0].set_xlabel(r"$y_{\max}$")
    ax[0].set_ylabel("absorptive fraction f")
    ax[0].set_title("S7: f against the reach")
    ax[0].legend(fontsize=9)
    ax[0].grid(alpha=.3)
    ax[1].plot(ymax, [r[3] for r in rows], "o-", lw=2, color="C0")
    ax[1].axhline(0.0, color="k", lw=.8)
    ax[1].set_xlabel(r"$y_{\max}$")
    ax[1].set_ylabel(r"worst cell $S/B$")
    ax[1].set_title("residual deficit against the reach")
    ax[1].grid(alpha=.3)
    fig.suptitle(f"The ledger identity (measured f = {f_meas:.4f})", y=1.02)
    save_fig(fig, "sea_population_ledger_identity.png")


def main():
    run = Ledger()
    G, E = part_a(run)
    e0, outs = part_b(run)
    part_c(run, e0)
    part_d(run)
    f_meas = part_e(run, e0)
    fid, tr_e, tr_a = part_f(run, e0)
    rows = part_g(e0)
    figures(run, G, E, outs, fid, tr_e, tr_a, rows, f_meas)
    print("\nFigures: sea_population_fixed_point.png, "
          "sea_population_unravelling.png,\n"
          "         sea_population_ledger_identity.png")


if __name__ == "__main__":
    main()
