"""Contact-generated kernel weights in the ledger -- companion to step 22 section 2.

The weight of channel q at x is an integral over leg half-separations y,

    K_q(x) = sum_j C[x, j, q],   C = Re[(1/n_p) M_res(x, s_j) e^{2 pi i j q/n_p}]

(the discrete form of -(B dp/hbar) int dy U_res w sin(2 xi_q y/hbar)).
Read as contacts: each sampled separation y_j is one encounter, carrying
the signed weight n_p C[x, j, q] for EVERY channel at once.

A  identity check: sum_j C = K, and the sea-weighted form with S = B.
B  test 1, contact noise.  Pinned sea.  Each cell draws Poisson contacts at
   rate R_c per unit time; the step's kernel is the unbiased estimate.
   Correlated worst case: all bodies in a cell share the contacts.
C  test 2, the measure.  No sampling noise.  Contacts drawn from the
   ACTUAL ledger sea: weight C[x, j, q] * S(x + y_j, p) / B, so the
   kernel depends on the parent's row as well.  The force stays with
   streaming (pinned B), so any first moment the reweighting creates is
   a spurious force.
Reference throughout: the exact mesh (Ledger.qle_step) on the same grid.

Usage::

    PYTHONPATH=src python3 -u src/demo_contact_kernel.py [RATE ...] [--comp r|min]
        [--sea S|blind]

--comp chooses the per-contact force removal of part B: along the force
kernel r, or the least-change projection along xi_q (min).

--sea blind moves the ledger sea under postulate (S') of step 23
(force_blind_sea.md): inertial advection, no drift in p.
"""
import argparse
import numpy as np
from demo_emission_and_absorption import Ledger, packet, B, MU

ap = argparse.ArgumentParser(description="step 22 contact-kernel probe")
ap.add_argument("rates", type=float, nargs="*", help="contact rates R_c")
ap.add_argument("--comp", choices=("r", "min"), default="min")
ap.add_argument("--sea", choices=("S", "blind"), default="S",
                help="sea transport: postulate (S), or (S') of step 23")
ARGS = ap.parse_args()
rng = np.random.default_rng(12345)


def banner(t):
    print("\n" + "=" * 72 + "\n" + t + "\n" + "=" * 72, flush=True)


run = Ledger(v0=1.0, a=1.0, n_r=192, r_half=24.0, n_p=64, dp=0.25)
run.sea_force = ARGS.sea == "S"      # (S) as published, or (S') of step 23
n_p, n_r = run.n_p, len(run.r)
Q = np.arange(1, n_p // 2)                         # channels, as in channels()
m_res = np.fft.fft(run.k_full, axis=1)             # recover M_res(x, s_j)
J = np.arange(n_p)
phase = np.exp(2j * np.pi * J[:, None] * Q[None, :] / n_p)       # (j, q)
C = np.real(m_res[:, :, None] * phase[None, :, :]) / n_p         # (x, j, q)
shift = np.round(run.y / run.dr).astype(int)       # y_j in whole x-cells


def kernel_from_sea(sea):
    """k[x, p, q] = sum_j C[x, j, q] * S(x + y_j, p) / B."""
    w = sea / B
    ws = np.stack([np.roll(w, -sh, axis=0) for sh in shift])     # (j, x, p)
    return np.einsum('jxp,xjq->xpq', ws, C, optimize=True)


def channels_k(up, um, sea, dt, k):
    """Ledger.channels, unclamped, with a kernel k[x, p or 1, q].

    Every species mask is taken at the PARENT's row and rolled with the
    deposit.  The published version read the mask at the destination row
    p +- q; with a kernel whose sign varies with p -- the sea-weighted
    kernel of part C -- that assigned some deposits to the wrong species
    and drifted Sum E to 1.0010 (step 22 open item L-SP8, diagnosed in
    step 23, force_blind_sea.md section 6).  For a kernel whose sign
    depends on x only the two forms agree exactly.
    """
    R = lambda f, s: np.roll(f, s, axis=1)
    for iq, q in enumerate(Q):
        kq = k[:, :, iq]
        lam, sg = np.abs(kq), np.sign(kq)
        for parent, sp in ((up, 1.0), (um, -1.0)):
            D = lam * parent * dt
            t = np.broadcast_to(sg * sp, D.shape)
            pos, neg = t > 0, ~(t > 0)
            capA = np.clip(np.where(pos, R(um, -q), R(up, -q)), 0.0, None)
            capB = np.clip(np.where(pos, R(up, q), R(um, q)), 0.0, None)
            A = np.minimum(D, np.minimum(capA, capB))
            Em = D - A
            um -= R(np.where(pos, A, 0.0), q)
            up -= R(np.where(neg, A, 0.0), q)
            up -= R(np.where(pos, A, 0.0), -q)
            um -= R(np.where(neg, A, 0.0), -q)
            sea += A
            up += R(np.where(pos, Em, 0.0), q)
            um += R(np.where(neg, Em, 0.0), q)
            um += R(np.where(pos, Em, 0.0), -q)
            up += R(np.where(neg, Em, 0.0), -q)
            sea -= Em
    return up, um, sea

right = run.r[:, None] > 0.0
T = lambda w: float(w[np.broadcast_to(right, w.shape)].sum() / w.sum())
p0, dt = 1.2, 0.02
e0 = packet(run, r0=-8.0, p0=p0, sr=2.0, sp=0.25)
nstep = int(round(16.0 / (p0 / MU) / dt))

# ---------------------------------------------------------------- A
banner("A  identity: the contact sum reproduces the kernel")
k_sum = C.sum(axis=1)
print(f"   max |sum_j C - K| / max|K| = "
      f"{np.abs(k_sum - run.k[:, Q]).max() / np.abs(run.k).max():.2e}")
k_b = kernel_from_sea(np.full((n_r, n_p), B))
print(f"   max |K(S = B) - K| / max|K| = "
      f"{np.abs(k_b - run.k[:, None, Q]).max() / np.abs(run.k).max():.2e}")
print(f"   shifts y_j/dr: min {shift.min()}, max {shift.max()}"
      f"  (reach {run.y_max:.3f}, dr {run.dr:.4f})")

# mesh reference, once
emesh = e0.copy()
for _ in range(nstep):
    emesh = run.qle_step(emesh, dt)
print(f"   mesh reference: T_E = {T(emesh):.4f}   ({nstep} steps, dt = {dt})")


def run_ledger(kernel_fn, label):
    up, um, sea = run.prepare(e0, 6.0)
    fm_max, dev_max = 0.0, 0.0
    for _ in range(nstep):
        up, um, sea = run.stream3(up, um, sea, .5 * dt)
        k = kernel_fn(sea)
        fm = np.abs((Q * run.dp * k).sum(axis=-1)).max()
        fm_max = max(fm_max, fm)
        dev_max = max(dev_max, np.abs(sea / B - 1.0).max())
        up, um, sea = channels_k(up, um, sea, dt, k)
        up, um, sea = run.stream3(up, um, sea, .5 * dt)
    E = up - um
    err = np.linalg.norm(E - emesh) / np.linalg.norm(emesh)
    print(f"   {label:26s} T_E {T(E):7.4f}   |E-mesh|/|mesh| {err:6.3f}"
          f"   sum E {E.sum():7.4f}   max|S/B-1| {dev_max:6.3f}"
          f"   max|first moment| {fm_max:.2e}", flush=True)
    return E


# ---------------------------------------------------------------- B
banner("B  test 1: contact noise, pinned sea (per-cell shared contacts)")
gmax = run.gamma_tot.max()
print(f"   max Gamma(x) = {gmax:.3f} events per parent per unit time")
k_exact = run.k[:, None, Q]
run_ledger(lambda sea: k_exact, "exact kernel")
rates = ARGS.rates or [3.0, 30.0, 300.0, 3000.0]


def noisy(rc):
    lam = rc * dt / n_p                      # Poisson mean per (x, j)
    def fn(sea):
        N = rng.poisson(lam, size=(n_r, n_p)).astype(float)
        kh = np.einsum('xj,xjq->xq', N, C) * (n_p / (rc * dt))
        return kh[:, None, :]
    return fn


for rc in rates:
    if ARGS.comp == 'min':
        break
    run_ledger(noisy(rc), f"R_c = {rc:g} contacts/time")

# per-contact compensation: remove each contact's own first moment along
# the compensation direction r (the kernel of i*s, i.e. a pure force).
# Since sum_j m1_j = 0 already, the mean kernel is unchanged exactly.
s_sym = 2.0 * np.pi * np.fft.fftfreq(n_p, d=run.dp)
s_sym[n_p // 2] = 0.0
wy = np.cos(np.pi * run.y / (2.0 * run.y_max)) ** 2
r = np.real(np.fft.ifft(1j * s_sym * wy))[Q]                    # (q,)
xiQ = Q * run.dp
m1_ref = (xiQ * r).sum()
m1_j = np.einsum('xjq,q->xj', C, xiQ)                            # (x, j)
if ARGS.comp == 'min':      # least-change projection, along xi_q
    C_comp = C - (m1_j / (xiQ ** 2).sum())[:, :, None] * xiQ[None, None, :]
else:                                    # along the force kernel r
    C_comp = C - (m1_j / m1_ref)[:, :, None] * r[None, None, :]
print(f"   per-contact force before: max|m1_j| = {np.abs(m1_j).max():.3e};"
      f"  after: {np.abs(np.einsum('xjq,q->xj', C_comp, xiQ)).max():.1e}")
print(f"   mean kernel unchanged: max|sum_j C_comp - K|/max|K| = "
      f"{np.abs(C_comp.sum(axis=1) - run.k[:, Q]).max()/np.abs(run.k).max():.1e}")
wmax = n_p * np.abs(C_comp).sum(axis=2).max()
print(f"   largest single-contact total weight n_p*sum_q|C| = {wmax:.1f}"
      f"  (so R_c must exceed ~{wmax:.0f} for one contact to trigger at most one event)")


def noisy_comp(rc):
    lam = rc * dt / n_p
    def fn(sea):
        N = rng.poisson(lam, size=(n_r, n_p)).astype(float)
        kh = np.einsum('xj,xjq->xq', N, C_comp) * (n_p / (rc * dt))
        return kh[:, None, :]
    return fn


for rc in rates:
    run_ledger(noisy_comp(rc), f"R_c = {rc:g}, compensated")

# ---------------------------------------------------------------- C
banner("C  test 2: contacts drawn from the actual ledger sea")
run_ledger(kernel_from_sea, "sea-weighted kernel")
