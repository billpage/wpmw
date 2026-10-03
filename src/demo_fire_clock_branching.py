#!/usr/bin/env python3
r"""
Deterministic integrate-and-fire clocks for Cyganski's pair branching
(docs/supplement/poisson_kicks_and_pair_branching.md, Parts M-N;
step 23's Proposition Q11 in docs/analysis/force_blind_sea.md).

Same problem as src/demo_event_driven_branching.py: V = omega^2 x^2/2 +
V1 sin(k x), hbar = m = 1, omega = 0.4, V1 = 0.5, k = 1, a Gaussian at rest at
x0 = 6 with sigma_x = 0.7; the quadratic part is an exact flow and the sine acts
through signed pair branching (Lemma T1).  Only the clock changes:

  thin      Poisson by thinning (exact; the supplement's reference clock)
  fire      every particle integrates its signed rate Gamma(x) = V1 cos(k x)
            along its flight and fires at each integer crossing, with the
            crossing's direction as the event's orientation (Q11); only the
            initial phases are random
  fire-abs  the same with |Gamma| (a gross clock)

Parts
  M  No annihilation, t = 4.5: rms L2 error of rho(x) against |psi|^2 and
     e * sqrt(N0) over N0, for the three clocks with the same initial samples.
     Constant e * sqrt(N0) means unbiased; its size is the noise.  Also the
     sub-step convergence of the fire clock.
  N  With cell annihilation (h_x = 0.1, h_p = 0.125), t = 9: the clock against
     the annihilation bias of Proposition X5.

Usage:  PYTHONPATH=src python3 -u src/demo_fire_clock_branching.py
        (WPMW_QUICK=1 for a smoke test)
"""

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wpmwlib.event_branching import EventBranching  # noqa: E402

OM, V1, K = 0.4, 0.5, 1.0
X0, SX = 6.0, 0.7
SP = 1.0 / (2 * SX)
QUICK = os.environ.get("WPMW_QUICK", "") == "1"


def banner(s):
    print("\n" + "=" * 72 + "\n" + s + "\n" + "=" * 72, flush=True)


def V(x):
    return 0.5 * OM ** 2 * x ** 2 + V1 * np.sin(K * x)


# reference: split-operator Schroedinger, as Part F of the event-driven demo
NS, LS = 4096, 60.0
xs = (np.arange(NS) - NS // 2) * (LS / NS)
dxs = LS / NS
ks = 2 * np.pi * np.fft.fftfreq(NS, d=dxs)
psi = ((2 * np.pi * SX ** 2) ** -0.25 * np.exp(-(xs - X0) ** 2 / (4 * SX ** 2))).astype(complex)
dts = 0.0025
hv, kin = np.exp(-0.5j * V(xs) * dts), np.exp(-0.5j * ks ** 2 * dts)
rho_s = {}
for i in range(1, int(round(9.0 / dts)) + 1):
    psi = hv * np.fft.ifft(kin * np.fft.fft(hv * psi))
    for t in (4.5, 9.0):
        if i == int(round(t / dts)):
            rho_s[t] = np.abs(psi) ** 2

edges = np.arange(-12.0, 12.0 + 1e-9, 0.1)
HX = 0.1


def ref_bins(rho):
    cum = np.concatenate([[0.0], np.cumsum(rho) * dxs])
    xe = np.concatenate([xs - dxs / 2, [xs[-1] + dxs / 2]])
    return np.diff(np.interp(edges, xe, cum)) / HX


REF = {t: ref_bins(rho_s[t]) for t in rho_s}


def mc(N0, T, seed, clock, annihilate=False, cell=(0.1, 0.125), dt_fire=0.005):
    rng = np.random.default_rng(seed)
    x = rng.normal(X0, SX, N0)
    p = rng.normal(0.0, SP, N0)
    s = np.ones(N0, np.int8)
    eb = EventBranching(OM, V1, K, rule="pair", clock=clock, cell=cell,
                        rng=np.random.default_rng(10_000 + seed), dt_fire=dt_fire)
    if annihilate:
        out, _ = eb.run(x, p, s, [T], 0.1)
        x, p, s = out[T]
    else:
        x, p, s = eb.advance(x, p, s, 0.0, T)
    h, _ = np.histogram(x, bins=edges, weights=s.astype(float))
    return h / N0 / HX, x.size / N0, eb.n_events / N0


def l2(r, T):
    return float(np.sqrt(np.sum((r - REF[T]) ** 2) * HX))


CLOCKS = ("thin", "fire", "fire-abs")

banner("PART M  --  no annihilation, t = 4.5: the clock's bias and noise")
print("  rms over 4 seeds; each seed draws the same initial ensemble for every clock")
NM = [6250, 25000, 100000] if not QUICK else [6250]
for N0 in NM:
    for clk in CLOCKS:
        t0 = time.time()
        es, pops, evs, rs = [], [], [], []
        for rep in range(4):
            r, pop, ev = mc(N0, 4.5, 100 + rep, clk)
            es.append(l2(r, 4.5)); pops.append(pop); evs.append(ev); rs.append(r)
        e = float(np.sqrt(np.mean(np.square(es))))
        em = l2(np.mean(rs, axis=0), 4.5)
        print(f"  N0 = {N0:7d}  {clk:8s}:  rms L2 = {e:.4f}   e*sqrt(N0) = {e * np.sqrt(N0):6.2f}"
              f"   L2 of the 4-seed mean = {em:.4f}   events / N0 = {np.mean(evs):6.2f}"
              f"   population / N0 = {np.mean(pops):5.1f}   ({time.time() - t0:.0f} s)", flush=True)
print("  sub-step convergence of the fire clock (N0 = 25000, same 4 seeds):")
for dtf in (0.02, 0.01, 0.005):
    es = [l2(mc(25000, 4.5, 100 + rep, "fire", dt_fire=dtf)[0], 4.5) for rep in range(4)]
    print(f"    dt_fire = {dtf:5.3f}:  rms L2 = {np.sqrt(np.mean(np.square(es))):.4f}", flush=True)

banner("PART N  --  with cell annihilation (0.1 x 0.125), t = 9")
NN = [(100000, 4), (400000, 2)] if not QUICK else [(25000, 2)]
for N0, R in NN:
    for clk in ("thin", "fire"):
        t0 = time.time()
        res = [mc(N0, 9.0, 1000 * rep + N0 % 997, clk, annihilate=True) for rep in range(R)]
        e = float(np.sqrt(np.mean([l2(r, 9.0) ** 2 for r, _, _ in res])))
        print(f"  N0 = {N0:7d}  {clk:5s}:  rms L2 = {e:.4f}   events / N0 = "
              f"{np.mean([ev for _, _, ev in res]):6.2f}   ({time.time() - t0:.0f} s)", flush=True)
