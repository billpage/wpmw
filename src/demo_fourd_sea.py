"""A four-dimensional sea ledger -- step 24, open item V-SP1.

Companion to ``docs/analysis/sea_depletion.md``.  Proposition V5 settled
two bodies with a pair potential from a product of centre-of-mass and
relative states: the sea pays exactly the relative motion's negativity.
This demo runs the case V5 does not cover: two particles on a line in a
non-quadratic external potential, interacting through a non-quadratic pair
potential, so that the centre of mass and the relative motion couple.

    H = p1^2/2 + p2^2/2 + V(x1) + V(x2) + Wp(x1 - x2)
    V(x)  = -V0 sech^2(x/a)            (step 24's Poeschl-Teller well)
    Wp(r) = w0 exp(-r^2 / 2 s^2)       (a soft repulsion)

The ledger is step 24's minimal ledger in four dimensions.  E(x1, x2, p1,
p2) is transported spectrally (x_i at p_i, p_i under the compensated
force); the sea S has density B^2 per unit of four-dimensional phase
space.  By Proposition M1 of step 21 the residual kernel of a sum of
one-body and pair terms has three families of channels:
    family 1   jumps p1 -> p1 +- xi, rate from V at x1;
    family 2   jumps p2 -> p2 +- xi, rate from V at x2;
    family 3   jumps (p1, p2) -> (p1 +- xi, p2 -+ xi), rate from Wp at
               r = x1 - x2 (relative momentum only).
Each family's kernel is the one-dimensional compensated kernel of step 16
(``SeaLedger``), built for V on the x grid and for Wp on the grid of
differences.  Events are allocated absorptive-first exactly as in one
dimension, a contact sink keeps u+ = E+ and u- = E- (Proposition Q8), and
rule (F) caps an emission at the sea in its parent's cell.

Checks printed: the QLE on the same mesh (spectral transport, exact
exponential kick by the residual symbols) against the ledger's E;
Proposition V2 in four dimensions (D = dM- - dN_tr/2); the floor's
threshold bit for bit (Proposition V3); and, for orientation, the exact
Schroedinger evolution on the (x1, x2) grid, with its entanglement (the
purity of one particle's reduced state) against the purity read from the
mesh Wigner function, Tr rho1^2 = h * integral W1^2.

Run as::

    WPMW_OUTPUT=... PYTHONPATH=src python3 -u src/demo_fourd_sea.py [--gpu]

``--gpu`` uses CuPy (a Kaggle GPU kernel); the default grid then has
48^4 cells.  ``--small`` uses 24^4 cells, for testing on a CPU.
Writes ``fourd_sea<tag>.csv`` (one row per run) and
``fourd_sea_trace<tag>.csv`` (lambda*(t) and the ledger identity) through
``output_path``.
"""

from __future__ import annotations

import argparse
import csv
import time

import numpy as np

import demo_sea_depletion as sd
from wpmwlib.wpmw_utils import output_path

HBAR = sd.HBAR
B = sd.B
B4 = B * B


def well(v0, a):
    return lambda r: -v0 / np.cosh(r / a) ** 2


def gauss(w0, s):
    return lambda r: w0 * np.exp(-r ** 2 / (2.0 * s * s))


# ----------------------------------------------------------------------
class FourD:
    def __init__(self, xp, n_x=48, dx=0.3125, n_p=48, dp=0.25,
                 v0=4.0, a=2.0, w0=2.0, s=1.0):
        self.xp = xp
        self.n_x, self.dx, self.n_p, self.dp = n_x, dx, n_p, dp
        L = n_x * dx
        ext = sd.SeaLedger(well(v0, a), n_r=n_x, r_half=0.5 * L, n_p=n_p,
                           dp=dp)
        rel = sd.SeaLedger(gauss(w0, s), n_r=2 * n_x, r_half=L, n_p=n_p,
                           dp=dp)
        self.x = ext.r                                   # x_i = -L/2 + i dx
        self.p = ext.p                                   # FFT order
        self.Vx = well(v0, a)(self.x)
        self.Wf = gauss(w0, s)
        i1, i2 = np.meshgrid(np.arange(n_x), np.arange(n_x), indexing="ij")
        J = i1 - i2 + n_x                                # r = (i1 - i2) dx
        assert np.allclose(rel.r[J], self.x[i1] - self.x[i2])
        self.y_max = ext.y_max
        self.vol = (dx * dp) ** 2
        Q = np.arange(1, n_p // 2)
        self.Q = Q
        to = xp.asarray
        # family kernels, signed, q >= 1
        self.k1 = to(ext.k[:, Q])                        # (n_x, nq) at x1
        self.kW = to(rel.k[J][:, :, Q])                  # (n_x, n_x, nq)
        dv_ext = ext.dv_eff
        dv_W = rel.dv_eff[J]
        self.dv1 = to(dv_ext[:, None] + dv_W)            # dp1/dt = -dv1
        self.dv2 = to(dv_ext[None, :] - dv_W)            # dp2/dt = -dv2
        kx = 2.0 * np.pi * np.fft.fftfreq(n_x, d=dx)
        sgrid = 2.0 * np.pi * np.fft.fftfreq(n_p, d=dp)
        P = self.p
        self.ph1 = to(kx[:, None, None, None] * P[None, None, :, None])
        self.ph2 = to(kx[None, :, None, None] * P[None, None, None, :])
        self.kx, self.s = kx, sgrid
        # residual symbol of the full jump operator, for the QLE reference
        k1i = np.arange(n_p)
        dif = (k1i[:, None] - k1i[None, :]) % n_p        # (s1 - s2) index
        sym1 = ext.sym_e                                 # (n_x, n_p)
        symW = rel.sym_e[J]                              # (n_x, n_x, n_p)
        M = (sym1[:, None, :, None] + sym1[None, :, None, :]
             + symW[:, :, dif])
        self.M = to(M)
        # observables
        X1 = np.broadcast_to(self.x[:, None, None, None],
                             (n_x, n_x, n_p, n_p))
        X2 = np.broadcast_to(self.x[None, :, None, None], X1.shape)
        P1 = np.broadcast_to(P[None, None, :, None], X1.shape)
        P2 = np.broadcast_to(P[None, None, None, :], X1.shape)
        Ham = (0.5 * P1 ** 2 + 0.5 * P2 ** 2
               + self.Vx[:, None, None, None] + self.Vx[None, :, None, None]
               + self.Wf(X1 - X2))
        self.obs = {k: to(v) for k, v in dict(x1=X1, x2=X2, p1=P1, p2=P2,
                                               x1x2=X1 * X2, H=Ham).items()}

    # -- arrays -----------------------------------------------------------
    def product_state(self, xa, xb, sig):
        sp = HBAR / (2.0 * sig)

        def g(x0):
            w = np.exp(-(self.x[:, None] - x0) ** 2 / (2 * sig ** 2)
                       - self.p[None, :] ** 2 / (2 * sp ** 2))
            return w / (w.sum() * self.dx * self.dp)

        w1, w2 = g(xa), g(xb)
        return self.xp.asarray(w1[:, None, :, None] * w2[None, :, None, :])

    def l1(self, E):
        return float(self.xp.abs(E).sum()) * self.vol

    def neg(self, E):
        return float(self.xp.maximum(-E, 0.0).sum()) * self.vol

    def moments(self, E):
        return {k: float((E * v).sum()) * self.vol
                for k, v in self.obs.items()}

    def purity1(self, E):
        """Tr rho1^2 = h * integral W1^2, W1 the reduced Wigner function."""
        W1 = E.sum(axis=(1, 3)) * self.dx * self.dp
        return float((W1 ** 2).sum()) * self.dx * self.dp * 2 * np.pi * HBAR

    # -- transport ----------------------------------------------------------
    def stream(self, F, h, feel=True):
        xp = self.xp
        F = xp.fft.ifft(xp.fft.fft(F, axis=0) * xp.exp(-1j * h * self.ph1),
                        axis=0)
        F = xp.fft.ifft(xp.fft.fft(F, axis=1) * xp.exp(-1j * h * self.ph2),
                        axis=1).real
        if feel:
            s = xp.asarray(self.s)
            F = xp.fft.ifft(xp.fft.fft(F, axis=2) * xp.exp(
                1j * h * self.dv1[:, :, None, None]
                * s[None, None, :, None]), axis=2)
            F = xp.fft.ifft(xp.fft.fft(F, axis=3) * xp.exp(
                1j * h * self.dv2[:, :, None, None]
                * s[None, None, None, :]), axis=3).real
        return F

    def qle_step(self, E, dt):
        xp = self.xp
        E = self.stream(E, 0.5 * dt)
        E = xp.fft.ifft2(xp.fft.fft2(E, axes=(2, 3)) * xp.exp(dt * self.M),
                         axes=(2, 3)).real
        return self.stream(E, 0.5 * dt)

    # -- events ---------------------------------------------------------------
    def channels(self, up, um, S, dt, rule, beta, absorb, margin):
        """The 1D allocation of ``SeaLedger.channels_supply``, for the three
        families.  A daughter offset is a tuple of (axis, shift)."""
        xp = self.xp
        n_abs = n_emi = blocked = 0.0
        for fam in (1, 2, 3):
            for j, q in enumerate(self.Q):
                if fam == 1:
                    kq = self.k1[:, j][:, None, None, None]
                    plus = ((2, q),)
                elif fam == 2:
                    kq = self.k1[:, j][None, :, None, None]
                    plus = ((3, q),)
                else:
                    kq = self.kW[:, :, j][:, :, None, None]
                    plus = ((2, q), (3, -q))
                lam = xp.abs(kq)
                if float(lam.max()) < 1e-14:
                    continue
                sg = xp.sign(kq)
                minus = tuple((ax, -sh) for ax, sh in plus)

                def at(F, off):                  # F at parent + off
                    for ax, sh in off:
                        F = xp.roll(F, -sh, axis=ax)
                    return F

                def put(F, off):                 # move F from parent to +off
                    for ax, sh in off:
                        F = xp.roll(F, sh, axis=ax)
                    return F

                for parent, spc in ((up, 1.0), (um, -1.0)):
                    D = lam * parent * dt
                    if float(D.max()) <= 0.0:
                        continue
                    t = sg * spc
                    pos = t > 0
                    if absorb:
                        capA = xp.where(pos, at(um, plus), at(up, plus))
                        capB = xp.where(pos, at(up, minus), at(um, minus))
                        A = xp.minimum(D, xp.minimum(capA, capB))
                    else:
                        A = xp.zeros_like(D)
                    Em = D - A
                    due = Em > 0
                    if margin is not None and bool(due.any()):
                        m = float(xp.where(due, S - Em, xp.inf).min()) / B4
                        margin[0] = min(margin[0], m)
                    if rule == "floor":
                        Er = xp.minimum(Em, xp.maximum(S, 0.0))
                        blocked += float((Em - Er).sum())
                        Em = Er
                    n_abs += float(A.sum())
                    n_emi += float(Em.sum())
                    aP, aM = put(A, plus), put(A, minus)
                    um -= xp.where(put(pos, plus), aP, 0.0)
                    up -= xp.where(put(pos, plus), 0.0, aP)
                    up -= xp.where(put(pos, minus), aM, 0.0)
                    um -= xp.where(put(pos, minus), 0.0, aM)
                    S += A
                    eP, eM = put(Em, plus), put(Em, minus)
                    up += xp.where(put(pos, plus), eP, 0.0)
                    um += xp.where(put(pos, plus), 0.0, eP)
                    um += xp.where(put(pos, minus), eM, 0.0)
                    up += xp.where(put(pos, minus), 0.0, eM)
                    S -= Em
                    xp.maximum(up, 0.0, out=up)
                    xp.maximum(um, 0.0, out=um)
        return up, um, S, n_abs, n_emi, blocked

    # -- one run ----------------------------------------------------------------
    def run(self, E0, t_max, dt, sea_force=True, beta=1.0, rule="none",
            absorb=True, ref=True, every=1.0, motion=None):
        """``motion``: a ``Collisions4D`` applied to the sea once per step
        after the events (the amended (S'); Part G of the 1D demo)."""
        xp = self.xp
        E = E0.copy()
        S = xp.full(E0.shape, beta * B4)
        S0 = float(S.sum()) * self.vol
        M0 = self.neg(E0)
        n_tr = 0.0
        rE = E0.copy() if ref else None
        margin = [beta]
        tot = np.zeros(4)
        tr = []
        n_steps = int(round(t_max / dt))
        k_every = max(1, int(round(every / dt)))
        t_run = time.time()
        for step in range(n_steps + 1):
            if step % k_every == 0 or step == n_steps:
                row = dict(t=step * dt, D=S0 - float(S.sum()) * self.vol,
                           dM=self.neg(E) - M0, ntr2=0.5 * n_tr,
                           lam=beta - margin[0],
                           worst=beta - float(S.min()) / B4,
                           purity=self.purity1(E),
                           pairs=float(S.sum()) * self.vol)
                if ref:
                    row["fid"] = float(xp.linalg.norm(E - rE)
                                       / xp.linalg.norm(rE))
                    row["purity_qle"] = self.purity1(rE)
                tr.append(row)
            if step == n_steps:
                break
            for half in (0, 1):
                if half:
                    up, um = xp.maximum(E, 0.0), xp.maximum(-E, 0.0)
                    up, um, S, x1, x2, x3 = self.channels(
                        up, um, S, dt, rule, beta, absorb, margin)
                    c = xp.minimum(up, um)
                    S += c
                    E = up - um
                    tot += (x1, x2, x3, float(c.sum()))
                    if motion is not None:
                        S = motion.step(S, dt)
                n0 = self.l1(E)
                E = self.stream(E, 0.5 * dt)
                n_tr += self.l1(E) - n0
                S = self.stream(S, 0.5 * dt, feel=sea_force)
            if ref:
                rE = self.qle_step(rE, dt)
            if step == 9:
                print(f"      (first 10 steps: {time.time() - t_run:.1f} s)",
                      flush=True)
        n_ev = tot[0] + tot[1]
        return dict(E=E, S=S, ref=rE, trace=tr, lam=beta - margin[0],
                    f=tot[0] / max(n_ev, 1e-30), n_ev=n_ev * self.vol,
                    blocked=tot[2] * self.vol,
                    blocked_frac=tot[2] / max(tot[1] + tot[2], 1e-30),
                    D=S0 - float(S.sum()) * self.vol, dM=self.neg(E) - M0,
                    ntr2=0.5 * n_tr)


# ----------------------------------------------------------------------
class Collisions4D:
    """Definition (C) of step 24 in four dimensions: Focus/Defocus between
    aligned pairs at one joint position (x1, x2), along the momentum
    directions of the three event families (Proposition M1),
        family 1   v = (q, 0)     rate from V at x1,
        family 2   v = (0, q)     rate from V at x2,
        family 3   v = (q, -q)    rate from Wp at x1 - x2,
    with J_v(n) = r_v (s_n^2 - s_{n-v} s_{n+v}) and
    ds_n/dt = sum_v (-2 J_v(n) + J_v(n - v) + J_v(n + v)).  Each direction
    is a one-dimensional lattice of Proposition V10, so pairs and pair
    momentum are conserved and the collision entropy never increases.
    kind 'collide-k': every channel, r_v B^2 = scale |K_v| / 2, the family's
                      kernel at the cell (the 1D kernel rate, density B^2);
         'collide-u': v = (1, 0) and (0, 1) only, r B^2 = scale.
    A cell in debt (s < 0, the sea without a floor) holds no pairs to
    collide: the products use max(s, 0), which changes nothing while s > 0.
    Heun's method, with substeps for the linear symbol's largest rate."""

    def __init__(self, fd, kind, scale=1.0):
        xp = fd.xp
        self.xp, self.kind = xp, kind
        self.chan = []                      # (shift tuple, rate array)
        if kind == "collide-k":
            for j, q in enumerate(fd.Q):
                q = int(q)
                r1 = scale * xp.abs(fd.k1[:, j]) / (2.0 * B4)
                if float(r1.max()) > 1e-14 / B4:
                    self.chan.append((((2, q),), r1[:, None, None, None]))
                    self.chan.append((((3, q),), r1[None, :, None, None]))
                rW = scale * xp.abs(fd.kW[:, :, j]) / (2.0 * B4)
                if float(rW.max()) > 1e-14 / B4:
                    self.chan.append((((2, q), (3, -q)),
                                      rW[:, :, None, None]))
        elif kind == "collide-u":
            r = xp.asarray(scale / B4)
            self.chan = [(((2, 1),), r), (((3, 1),), r)]
        else:
            raise ValueError(kind)
        self.lam = float(sum(16.0 * float(xp.max(r)) * B4
                             for _, r in self.chan))

    def rhs(self, S):
        xp = self.xp
        Sp = xp.maximum(S, 0.0)
        out = xp.zeros_like(S)
        for sh, r in self.chan:
            ax = tuple(a for a, _ in sh)
            fwd = tuple(k for _, k in sh)
            bwd = tuple(-k for _, k in sh)
            # S at n - v and n + v
            Sm = xp.roll(Sp, fwd, axis=ax)
            Sq = xp.roll(Sp, bwd, axis=ax)
            J = r * (Sp * Sp - Sm * Sq)
            out += -2.0 * J + xp.roll(J, fwd, axis=ax) + xp.roll(J, bwd,
                                                                 axis=ax)
        return out

    def step(self, S, h):
        n = max(1, int(np.ceil(h * self.lam / 0.5)))
        hs = h / n
        for _ in range(n):
            k1 = self.rhs(S)
            k2 = self.rhs(S + hs * k1)
            S = S + 0.5 * hs * (k1 + k2)
        return S


# ----------------------------------------------------------------------
def schroedinger(fd, xa, xb, sig, t_max, dt, every=1.0):
    """Split-operator on the (x1, x2) grid of the mesh, no window: the
    exact evolution the mesh QLE approximates.  Returns moments and the
    purity of particle 1's reduced state at the trace times."""
    n, dx = fd.n_x, fd.dx
    x = fd.x
    g = lambda x0: np.exp(-(x - x0) ** 2 / (4 * sig ** 2))
    psi = (g(xa)[:, None] * g(xb)[None, :]).astype(complex)
    psi /= np.sqrt((np.abs(psi) ** 2).sum() * dx * dx)
    k = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    K = 0.5 * (k[:, None] ** 2 + k[None, :] ** 2)
    U = (fd.Vx[:, None] + fd.Vx[None, :]
         + fd.Wf(x[:, None] - x[None, :]))
    ek = np.exp(-1j * K * dt / HBAR)
    eu = np.exp(-0.5j * U * dt / HBAR)
    out = []
    n_steps = int(round(t_max / dt))
    k_every = max(1, int(round(every / dt)))
    for step in range(n_steps + 1):
        if step % k_every == 0 or step == n_steps:
            rho = np.abs(psi) ** 2 * dx * dx
            sv = np.linalg.svd(psi * dx, compute_uv=False)
            out.append(dict(t=step * dt,
                            x1=float((rho * x[:, None]).sum()),
                            x2=float((rho * x[None, :]).sum()),
                            x1x2=float((rho * x[:, None] * x[None, :]).sum()),
                            purity=float((sv ** 4).sum())))
        if step == n_steps:
            break
        psi = eu * psi
        psi = np.fft.ifft2(ek * np.fft.fft2(psi))
        psi = eu * psi
    return out


# ----------------------------------------------------------------------
# lambda*(T = 12) of the one-dimensional Poeschl-Teller well, step 24 section
# 4 (demo_sea_depletion.py, Part C), for comparison
WELL_1D = {"Sabs": 1.295, "Spabs": 0.765, "Spemi": 0.994}
MODE_STYLE = {"Sabs": ("#2a78d6", "-", "(S), absorptive first"),
              "Spabs": ("#eb6834", "--", "(S′), absorptive first"),
              "Spemi": ("#1baf7a", "-.", "(S′), emissive")}


def plot(trace_csv, sch_csv, fig_name="fourd_sea.png"):
    """Figure from a finished run's CSVs: lambda*(t) per mode against the
    one-dimensional well, and the purity of one particle's state from the
    ledger, the mesh QLE and Schroedinger."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from wpmwlib.wpmw_utils import docs_path

    tr = {}
    with open(trace_csv) as fh:
        for r in csv.DictReader(fh):
            tr.setdefault(r["mode"], []).append(
                {k: float(v) for k, v in r.items() if k != "mode"})
    with open(sch_csv) as fh:
        sch = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(fh)]
    ts = np.array([r["t"] for r in sch])
    ps = np.array([r["purity1"] for r in sch])
    first = next(iter(tr.values()))
    pq = np.array([r["purity1_qle"] for r in first])
    tq = np.array([r["t"] for r in first])
    bad = np.abs(pq - np.interp(tq, ts, ps)) > 0.05 * np.interp(tq, ts, ps)
    t_res = float(tq[np.argmax(bad)]) if bad.any() else None

    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.4))
    for m, rs in tr.items():
        col, ls, lab = MODE_STYLE[m]
        t = [r["t"] for r in rs]
        ax[0].plot(t, [r["lambda_star"] for r in rs], ls, color=col, lw=2,
                   label=f"4D, {lab}")
        ax[0].plot([t[-1] + 0.25], [WELL_1D[m]], "D", color=col, ms=8,
                   mfc="white", mew=2)
    ax[0].plot([], [], "D", color="0.35", mfc="white", mew=2,
               label="1D well at $t = 12$")
    ax[0].axhline(1.0, color="0.45", lw=1, ls=":")
    ax[0].text(8.4, 0.86, "physical density", color="0.3", fontsize=9)
    ax[0].set_xlabel("$t$")
    ax[0].set_ylabel(r"$\lambda^*(t)$, in units of $B^2$")
    ax[0].set_title("Depth the floor needs: two bodies against one")
    ax[0].legend(fontsize=9, loc="upper left")

    m0 = "Sabs" if "Sabs" in tr else next(iter(tr))
    ax[1].plot(ts, ps, color="0.1", lw=2, label="Schrödinger (exact)")
    ax[1].plot(tq, pq, "--", color="#eb6834", lw=2, label="mesh QLE")
    ax[1].plot([r["t"] for r in tr[m0]], [r["purity1"] for r in tr[m0]],
               "o", color="#2a78d6", ms=8, mfc="none", mew=2,
               label="minimal ledger $E$")
    ax[1].set_xlabel("$t$")
    ax[1].set_ylabel(r"purity Tr $\rho_1^2$")
    ax[1].set_title("Entanglement made by the pair potential")
    ax[1].legend(fontsize=9)
    for a_ in ax:
        a_.grid(alpha=0.25)
        if t_res is not None:
            a_.axvspan(t_res, ts[-1] + 0.5, color="0.5", alpha=0.08, lw=0)
    if t_res is not None:
        for a_, y in ((ax[0], 0.17), (ax[1], 0.93)):
            a_.text(t_res + 0.2, y, "mesh off Schrödinger\nby > 5 %",
                    color="0.3", fontsize=9, va="top",
                    transform=a_.get_xaxis_transform() if a_ is ax[0]
                    else a_.transData)
    fig.tight_layout()
    fig.savefig(output_path(fig_name), dpi=150, bbox_inches="tight")
    dp = docs_path(fig_name)
    if dp:
        fig.savefig(dp, dpi=150, bbox_inches="tight")
    print(f"wrote {output_path(fig_name)}; mesh QLE purity first off"
          f" Schroedinger by > 5 % at t = {t_res}")


def collisions_check(fd, xp):
    """Collisions alone on a ripple of the sea: pairs and pair momentum
    conserved, collision entropy decreasing (Proposition V10 in 4D)."""
    P = np.asarray(fd.p)
    shape = (fd.n_x, fd.n_x, fd.n_p, fd.n_p)
    rng = np.random.default_rng(3)
    S = np.full(shape, B4)
    core = (np.abs(P) < 1.5)                    # away from the momentum ends
    bump = rng.uniform(-0.5, 0.5, shape) * core[None, None, :, None] \
        * core[None, None, None, :]
    S = xp.asarray(S * (1.0 + bump))
    p1 = xp.asarray(np.broadcast_to(P[None, None, :, None], shape))
    p2 = xp.asarray(np.broadcast_to(P[None, None, None, :], shape))

    def stats(S):
        x = S / B4
        return (float(S.sum()), float((S * p1).sum()), float((S * p2).sum()),
                float((x * xp.log(x) - x + 1.0).sum()))
    for kind in ("collide-k", "collide-u"):
        mo = Collisions4D(fd, kind)
        if kind == "collide-k":
            # the momentum grid is periodic: a long jump from the ripple
            # would wrap and spoil the momentum sum, so the check keeps the
            # channels with q <= 2 (pairs and entropy hold for all)
            full = len(mo.chan)
            mo.chan = [c for c in mo.chan if max(abs(k) for _, k in c[0])
                       <= 2]
            print(f"  collide-k check: {len(mo.chan)} of {full} channels"
                  " (q <= 2)")
        S1 = S.copy()
        s0 = stats(S1)
        out = [s0]
        for _ in range(5):
            S1 = mo.step(S1, 0.02)
            out.append(stats(S1))
        o = np.array(out)
        print(f"  {kind}: {len(mo.chan)} channels; after 5 steps"
              f" d(pairs)/pairs {(o[-1, 0] - o[0, 0]) / o[0, 0]:+.1e},"
              f" d(sum p1 S)/pairs {(o[-1, 1] - o[0, 1]) / o[0, 0]:+.1e},"
              f" d(sum p2 S)/pairs {(o[-1, 2] - o[0, 2]) / o[0, 0]:+.1e},"
              f" entropy {o[0, 3]:.4g} -> {o[-1, 3]:.4g}, decreasing at every"
              f" step: {bool(np.all(np.diff(o[:, 3]) < 0))}", flush=True)


def collisions_study(fd, E0, a, xp):
    """lambda*(t) of the two coupled bodies under the amended (S'), for
    each mode and sea motion; E checked against the run without collisions
    (Proposition V8)."""
    collisions_check(fd, xp)
    modes = {"Sabs": (True, True), "Spabs": (False, True),
             "Spemi": (False, False)}
    rows, trows = [], []

    def dump():
        with open(output_path(f"fourd_sea_collisions{a.tag}.csv"), "w",
                  newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["mode", "motion", "lambda_star", "f", "events",
                        "D", "dMneg", "max_dE_vs_none", "seconds"])
            w.writerows(rows)
        with open(output_path(f"fourd_sea_collisions_trace{a.tag}.csv"), "w",
                  newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["mode", "motion", "t", "lambda_star", "worst", "D",
                        "pairs", "purity1"])
            w.writerows(trows)

    for name in a.modes.split(","):
        force, absorb = modes[name]
        E_none = None
        for m in a.motions:
            kind, _, sc = m.partition(":")
            mo = None if kind == "none" else Collisions4D(fd, kind,
                                                         float(sc or 1.0))
            t1 = time.time()
            print(f"\n  {name}, sea motion {m}:"
                  + (f" {len(mo.chan)} channels, {mo.lam:.1f} max rate"
                     if mo else ""), flush=True)
            o = fd.run(E0, a.t_max, a.dt, sea_force=force, absorb=absorb,
                       ref=False, motion=mo)
            if kind == "none":
                E_none = o["E"]
            dE = (float(xp.abs(o["E"] - E_none).max()) if E_none is not None
                  else float("nan"))
            wall = time.time() - t1
            print(f"    lambda* = {o['lam']:.4f}, f = {o['f']:.4f},"
                  f" D = {o['D']:.5f}, max|E - E(none)| = {dE:.1e}"
                  f"  ({wall:.0f} s)")
            print("    t, lambda*(t), worst cell, D, sea pairs:")
            for r in o["trace"]:
                print(f"      {r['t']:5.1f} {r['lam']:7.4f} {r['worst']:7.4f}"
                      f" {r['D']:8.5f} {r['pairs']:10.5f}", flush=True)
                trows.append([name, m, r["t"], r["lam"], r["worst"], r["D"],
                              r["pairs"], r["purity"]])
            rows.append([name, m, o["lam"], o["f"], o["n_ev"], o["D"],
                         o["dM"], dE, wall])
            dump()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--gpu", action="store_true")
    ap.add_argument("--small", action="store_true",
                    help="24^4 cells (dx 0.625, dp 0.5), for testing")
    ap.add_argument("--t-max", type=float, default=12.0)
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--w0", type=float, default=2.0)
    ap.add_argument("--modes", default="Sabs,Spabs,Spemi")
    ap.add_argument("--betas", type=float, nargs="*", default=[0.5, 1.0],
                    help="floor depths run for every mode, besides the"
                         " threshold check at 0.98 and 1.02 lambda*")
    ap.add_argument("--tag", default="")
    ap.add_argument("--motions", nargs="+", metavar="KIND[:SCALE]",
                    help="amended (S'): run each mode under these sea"
                         " motions (none, collide-k, collide-u) instead of"
                         " the floor study; writes fourd_sea_collisions*.csv")
    ap.add_argument("--plot", nargs=2, metavar=("TRACE_CSV", "SCHR_CSV"),
                    help="only draw the figure from a finished run's CSVs")
    a = ap.parse_args()
    if a.plot:
        plot(*a.plot)
        return
    if a.gpu:
        try:
            import cupy as xp
        except ImportError:                 # not on every Kaggle image
            import subprocess
            import sys
            subprocess.run([sys.executable, "-m", "pip", "install", "-q",
                            "cupy-cuda12x"], check=True)
            import cupy as xp
        print(f"  GPU: {xp.cuda.runtime.getDeviceProperties(0)['name']}")
    else:
        xp = np
    grid = (dict(n_x=24, dx=0.625, n_p=24, dp=0.5) if a.small
            else dict(n_x=48, dx=0.3125, n_p=48, dp=0.25))
    t0 = time.time()
    fd = FourD(xp, w0=a.w0, **grid)
    xa, xb, sig = 2.0, -2.0, 0.6
    E0 = fd.product_state(xa, xb, sig)
    print(f"4D sea ledger: grid {grid}, y_max {fd.y_max:.3f}, T {a.t_max},"
          f" dt {a.dt}, w0 {a.w0}, backend {'cupy' if a.gpu else 'numpy'}")
    print(f"  norm {float(E0.sum()) * fd.vol:.12f},"
          f" max W0 / B^2 = {float(E0.max()) / B4:.4f},"
          f" purity {fd.purity1(E0):.6f}", flush=True)

    if a.motions:
        collisions_study(fd, E0, a, xp)
        print(f"\nDone in {time.time() - t0:.0f} s.")
        return
    sch = schroedinger(fd, xa, xb, sig, a.t_max, 0.005)
    print("\n  Schroedinger (exact on the x grid):  t, <x1>, <x2>, cov,"
          " purity1")
    for r in sch:
        cov = r["x1x2"] - r["x1"] * r["x2"]
        print(f"    {r['t']:5.1f} {r['x1']:+8.4f} {r['x2']:+8.4f}"
              f" {cov:+8.4f} {r['purity']:8.5f}")

    def dump():
        """Write all three CSVs; called after every run, so a kernel
        that hits its time limit still leaves what it finished."""
        with open(output_path(f"fourd_sea{a.tag}.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["mode", "beta", "lambda_star", "f", "events",
                        "blocked_frac", "eps", "D", "dMneg", "dNtr_half",
                        "fid_qle", "x1", "x2", "p1", "p2", "H", "dP", "dH",
                        "dMneg_vs_base"])
            w.writerows(rows)
        with open(output_path(f"fourd_sea_trace{a.tag}.csv"), "w",
                  newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["mode", "t", "lambda_star", "worst", "D", "dMneg",
                        "dNtr_half", "fid_qle", "purity1", "purity1_qle"])
            w.writerows(trows)
        with open(output_path(f"fourd_sea_schroedinger{a.tag}.csv"), "w",
                  newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["t", "x1", "x2", "x1x2", "purity1"])
            w.writerows([[r["t"], r["x1"], r["x2"], r["x1x2"], r["purity"]]
                         for r in sch])

    modes = {"Sabs": (True, True), "Spabs": (False, True),
             "Spemi": (False, False)}
    rows, trows = [], []
    for name in a.modes.split(","):
        force, absorb = modes[name]
        t1 = time.time()
        base = fd.run(E0, a.t_max, a.dt, sea_force=force, absorb=absorb,
                      ref=True)
        lam = base["lam"]
        mo = fd.moments(base["E"])
        mr = fd.moments(base["ref"])
        last = base["trace"][-1]
        print(f"\n  {name}: lambda* = {lam:.4f}, f = {base['f']:.4f},"
              f" events {base['n_ev']:.3f}, |E - QLE| = {last['fid']:.3e}"
              f"  ({time.time() - t1:.0f} s)")
        print(f"    V2: D = {base['D']:.5f}, dM- = {base['dM']:.5f},"
              f" dN_tr/2 = {base['ntr2']:+.5f}, residual"
              f" {base['D'] - base['dM'] + base['ntr2']:.1e}")
        print(f"    purity1: ledger {last['purity']:.5f},"
              f" mesh QLE {last['purity_qle']:.5f},"
              f" Schroedinger {sch[-1]['purity']:.5f};"
              f"  <x1> {mo['x1']:+.4f} (QLE {mr['x1']:+.4f},"
              f" exact {sch[-1]['x1']:+.4f})")
        print("    t, lambda*(t), worst cell, D, dM-, |E-QLE|, purity1:")
        for r in base["trace"]:
            print(f"      {r['t']:5.1f} {r['lam']:7.4f} {r['worst']:7.4f}"
                  f" {r['D']:8.5f} {r['dM']:8.5f} {r['fid']:9.2e}"
                  f" {r['purity']:8.5f}", flush=True)
            trows.append([name, r["t"], r["lam"], r["worst"], r["D"],
                          r["dM"], r["ntr2"], r["fid"], r["purity"],
                          r["purity_qle"]])
        rows.append([name, "inf", lam, base["f"], base["n_ev"], 0.0, 0.0,
                     base["D"], base["dM"], base["ntr2"], last["fid"],
                     *[mo[k] for k in ("x1", "x2", "p1", "p2", "H")], 0.0,
                     0.0, 0.0])
        dump()
        for beta in sorted(set(a.betas + [round(0.98 * lam, 4),
                                          round(1.02 * lam, 4)])):
            t1 = time.time()
            o = fd.run(E0, a.t_max, a.dt, sea_force=force, absorb=absorb,
                       beta=beta, rule="floor", ref=False)
            m = fd.moments(o["E"])
            eps = float(xp.linalg.norm(o["E"] - base["E"])
                        / xp.linalg.norm(base["E"]))
            d = {k: m[k] - mo[k] for k in m}
            print(f"    beta {beta:7.4f}: blocked {o['blocked_frac']:.2e},"
                  f" eps {eps:.3e}, d<p1+p2> {d['p1'] + d['p2']:+.1e},"
                  f" dH {d['H']:+.1e}, dM- {o['dM'] - base['dM']:+.1e}"
                  f"  ({time.time() - t1:.0f} s)", flush=True)
            rows.append([name, beta, lam, o["f"], o["n_ev"],
                         o["blocked_frac"], eps, o["D"], o["dM"], o["ntr2"],
                         "", *[m[k] for k in ("x1", "x2", "p1", "p2", "H")],
                         d["p1"] + d["p2"], d["H"], o["dM"] - base["dM"]])
            dump()
    dump()
    print(f"\nDone in {time.time() - t0:.0f} s.")


if __name__ == "__main__":
    main()
