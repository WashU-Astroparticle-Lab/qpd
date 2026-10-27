"""Blind LED phase: estimate ``first_pulse_time_s`` from the waveform when only the LED period is known.

Each LED pulse redraws the offset charge. In a correctly cut LED window the
trace therefore visits two fixed I/Q points, one per parity; a cut in the
wrong place mixes in a stretch of the neighbouring charge.

Both methods are waveform-only and start from the same coarse stage: the
phases whose cuts make every window's 2 ms-smoothed trace the cleanest
two-cluster cloud (two-means profile on a 1 ms grid).

* ``'hmm'`` (default) refines the best coarse candidates by an exact
  change-point profile of the pinned-rate parity x quiet/burst chain at every
  onset. One forward and one backward pass per onset give the likelihood of
  the switch from the previous window's centres to the next window's at every
  sample, and these profiles are summed over onsets. The candidate with the
  highest whole-record HMM log-likelihood wins.
* ``'v1'`` refines with a per-sample, memoryless mixture change-point
  likelihood instead. It is faster and less accurate.
"""
from math import log

import numpy as np

from qpd.mlebench.led_hmm import smoothed_fit
from qpd.reconstruction.burst_hmm import burst_transition_matrix


# --- coarse stage: two-cluster profile over candidate phases -------------------------------------------------------
def _two_means_cost(pts):
    n = pts.size
    if n < 4: return 0.0
    x, y = pts.real-pts.real.mean(), pts.imag-pts.imag.mean()
    sxx, syy, sxy = x@x, y@y, x@y; th = 0.5*np.arctan2(2*sxy, sxx-syy); c, s = np.cos(th), np.sin(th)
    a = np.sort(c*x+s*y); perp = float(np.sum((-s*x+c*y)**2))
    cs, cs2 = np.cumsum(a), np.cumsum(a*a); k = np.arange(1, n)
    return perp+float(np.min(cs2[:-1]-cs[:-1]**2/k+(cs2[-1]-cs2[:-1])-(cs[-1]-cs[:-1])**2/(n-k)))


def _two_means_centres(pts):
    x, y = pts.real-pts.real.mean(), pts.imag-pts.imag.mean()
    th = 0.5*np.arctan2(2*(x@y), x@x-y@y); proj = np.cos(th)*x+np.sin(th)*y; o = np.argsort(proj); a = proj[o]; n = a.size
    cs, cs2 = np.cumsum(a), np.cumsum(a*a); k = np.arange(1, n)
    j = int(np.argmin(cs2[:-1]-cs[:-1]**2/k+(cs2[-1]-cs2[:-1])-(cs[-1]-cs[:-1])**2/(n-k)))+1
    return pts[o[:j]].mean(), pts[o[j:]].mean()


def profile_cost(y, offs, phi, P, m):
    pos = (offs-phi) % P; keep = pos <= P-m
    win = np.floor((offs-phi)/P).astype(int); y, win = y[keep], win[keep]
    parts = np.split(y, np.flatnonzero(np.diff(win))+1)
    return sum(_two_means_cost(p) for p in parts)/max(y.size, 1)


def prepare(iq, m=200):
    ys = np.convolve(iq, np.ones(m)/m, mode='valid'); d = np.diff(iq)
    sig2 = (0.5*(np.median(np.abs(d.real))+np.median(np.abs(d.imag)))/(0.67449*np.sqrt(2)))**2
    return ys, np.arange(ys.size), sig2


def coarse_phase(ys, off, P, m, step):
    grid = np.arange(0, P, step); cost = np.array([profile_cost(ys[::step], off[::step], g, P, m) for g in grid])
    med = np.median(cost); depth = (med-cost.min())/max(np.median(np.abs(cost-med))*1.4826, 1e-30)
    return int(grid[np.argmin(cost)]) % P, float(depth)


# --- 'v1' method: memoryless mixture change-point likelihood ------------------------------------------------------
def fine_curve(iq, ys, off, phi, P, m, h, sig2, clip):
    n = iq.size; cuts = np.arange(phi, n, P); edges = np.r_[0, cuts, n]; cen = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (off >= lo+h) & (off+m <= hi-h); pts = ys[sel][::10]
        cen.append(_two_means_centres(pts) if pts.size >= 20 else None)
    curves = []
    for k, c in enumerate(cuts):
        pre, post = cen[k], cen[k+1]
        if pre is None or post is None or c-h-1 < 0 or c+h+1 > n: continue
        z = iq[c-h-1:c+h+1]
        lp = np.logaddexp(-np.abs(z-pre[0])**2/(2*sig2), -np.abs(z-pre[1])**2/(2*sig2))
        lq = np.logaddexp(-np.abs(z-post[0])**2/(2*sig2), -np.abs(z-post[1])**2/(2*sig2))
        j = np.arange(1, 2*h+2); curves.append(np.cumsum(lp)[j-1]+(np.sum(lq)-np.cumsum(lq)[j-1]))
    if not curves: return None
    C = np.array(curves)
    if clip is not None: C = np.maximum(C-C.max(axis=1, keepdims=True), -clip)
    return C.sum(0)


def estimate(iq, fs, period_s, coarse=None, coarse_step=100, n_iter=1, clip=None, m=200, half_s=3e-3, prep=None):
    """'v1' estimate: dict(fp = first onset in [0, period) [s], depth = coarse-profile depth, phi0 = coarse phase [samples],
    lo/hi1, 2, 4 = offsets [s] around the fine maximum where the summed likelihood is within 1, 2, 4 nats)."""
    P = int(round(period_s*fs)); h = int(round(half_s*fs)); ys, off, sig2 = prep or prepare(iq, m)
    phi, depth = coarse or coarse_phase(ys, off, P, m, coarse_step); phi0 = phi; total = None
    for _ in range(n_iter):
        t = fine_curve(iq, ys, off, phi, P, m, h, sig2, clip)
        if t is None or not np.any(t): break
        total = t; phi = (phi-h+int(np.argmax(t))) % P
    out = dict(fp=((phi-0.5)/fs) % (P/fs), depth=depth, phi0=phi0)
    for dl in (1, 2, 4):
        if total is None: out[f'lo{dl}'], out[f'hi{dl}'] = -h/fs, h/fs; continue
        idx = np.flatnonzero(total >= total.max()-dl); best = int(np.argmax(total))
        out[f'lo{dl}'], out[f'hi{dl}'] = (idx.min()-best)/fs, (idx.max()-best)/fs
    return out


# --- 'hmm' method: exact HMM change-point profile (default) ---------------------------------------------------------


def forward_loglik4(log_emit4, trans):
    """Total log-likelihood of a 4-state chain (scaled forward recursion, scalar arithmetic for speed)."""
    peak = log_emit4.max(axis=0); e = np.exp(log_emit4-peak); ls = float(peak.sum())
    e0, e1, e2, e3 = (r.tolist() for r in e); (t00, t01, t02, t03), (t10, t11, t12, t13), (t20, t21, t22, t23), (t30, t31, t32, t33) = trans.tolist()
    a0, a1, a2, a3 = .25*e0[0], .25*e1[0], .25*e2[0], .25*e3[0]; s = a0+a1+a2+a3; ls += log(s); a0, a1, a2, a3 = a0/s, a1/s, a2/s, a3/s
    for k in range(1, len(e0)):
        b0 = (a0*t00+a1*t10+a2*t20+a3*t30)*e0[k]; b1 = (a0*t01+a1*t11+a2*t21+a3*t31)*e1[k]
        b2 = (a0*t02+a1*t12+a2*t22+a3*t32)*e2[k]; b3 = (a0*t03+a1*t13+a2*t23+a3*t33)*e3[k]
        s = b0+b1+b2+b3; ls += log(s); a0, a1, a2, a3 = b0/s, b1/s, b2/s, b3/s
    return ls


def record_loglik(iq, fs, phi, P, m=200, quiet_hz=20., p_burst=.3, tau=3.7e-3):
    """Whole-record log-likelihood with cuts at phi + k P: per window, 2 ms smoothed-fit 2-D centres, robust noise and the
    pinned-rate 4-state chain; windows shorter than 2m are scored as one 2-D Gaussian. Every sample is scored once."""
    n = iq.size; dt = 1/fs; edges = np.unique(np.r_[0, np.arange(phi % P, n, P), n]); total = 0.
    trans = burst_transition_matrix(max(quiet_hz/fs, 1e-9), p_burst, dt*fs/P, dt/tau)
    for lo, hi in zip(edges[:-1], edges[1:]):
        z = iq[lo:hi]
        if z.size < 2: continue
        d = np.diff(z); s2 = (0.5*(np.median(np.abs(d.real))+np.median(np.abs(d.imag)))/(0.67449*np.sqrt(2)))**2; s2 = max(s2, 1e-30)
        if z.size < 2*m:
            total += float(np.sum(-np.abs(z-z.mean())**2/(2*s2)))-z.size*log(2*np.pi*s2); continue
        fit, _ = smoothed_fit(z, m); ca, cb = fit.origin+fit.direction*fit.mu_a, fit.origin+fit.direction*fit.mu_b
        le2 = np.stack([-np.abs(z-ca)**2/(2*s2), -np.abs(z-cb)**2/(2*s2)])-log(2*np.pi*s2)
        total += forward_loglik4(le2[[0, 1, 0, 1]], trans)
    return total


def _emissions(z, ca, cb, s2):
    le = np.stack([-np.abs(z-ca)**2/(2*s2), -np.abs(z-cb)**2/(2*s2)])[[0, 1, 0, 1]]
    peak = le.max(axis=0); return np.exp(le-peak), peak


def switch_profile(z, pre, post, s2, trans):
    """Exact log-likelihood of z (4-state chain) with pre-model emissions before sample j and post-model from j on, j = 1..L-1."""
    L = z.size; ep, pp = _emissions(z, *pre, s2); eq, pq = _emissions(z, *post, s2); T = trans.tolist()
    (t00, t01, t02, t03), (t10, t11, t12, t13), (t20, t21, t22, t23), (t30, t31, t32, t33) = T
    A = np.empty(L); al = np.empty((L, 4)); e0, e1, e2, e3 = (r.tolist() for r in ep)
    a = [.25*e0[0], .25*e1[0], .25*e2[0], .25*e3[0]]; s = sum(a); a = [x/s for x in a]; acc = log(s)+pp[0]; al[0] = a; A[0] = acc
    for t in range(1, L):
        a0, a1, a2, a3 = a
        b0 = (a0*t00+a1*t10+a2*t20+a3*t30)*e0[t]; b1 = (a0*t01+a1*t11+a2*t21+a3*t31)*e1[t]
        b2 = (a0*t02+a1*t12+a2*t22+a3*t32)*e2[t]; b3 = (a0*t03+a1*t13+a2*t23+a3*t33)*e3[t]
        s = b0+b1+b2+b3; a = [b0/s, b1/s, b2/s, b3/s]; acc += log(s)+pp[t]; al[t] = a; A[t] = acc
    B = np.empty(L); be = np.empty((L, 4)); f0, f1, f2, f3 = (r.tolist() for r in eq)
    b = [1., 1., 1., 1.]; acc = 0.; be[L-1] = b; B[L-1] = 0.
    for t in range(L-2, -1, -1):
        v0, v1, v2, v3 = f0[t+1]*b[0], f1[t+1]*b[1], f2[t+1]*b[2], f3[t+1]*b[3]
        c0 = t00*v0+t01*v1+t02*v2+t03*v3; c1 = t10*v0+t11*v1+t12*v2+t13*v3
        c2 = t20*v0+t21*v1+t22*v2+t23*v3; c3 = t30*v0+t31*v1+t32*v2+t33*v3
        s = c0+c1+c2+c3; b = [c0/s, c1/s, c2/s, c3/s]; acc += log(s)+pq[t+1]; be[t] = b; B[t] = acc
    q = al[:-1] @ np.asarray(T); r = eq.T[1:]*be[1:]
    return A[:-1]+pq[1:]+B[1:]+np.log(np.sum(q*r, axis=1))   # index j-1 for switch at j = 1..L-1


def hmm_fine(iq, fs, phi, P, W_s, m=200, quiet_hz=20., p_burst=.3, tau=3.7e-3):
    """Best phase near phi: summed per-onset switch profiles over offsets -W..+W (sample resolution)."""
    n = iq.size; W = int(round(W_s*fs)); dt = 1/fs; trans = burst_transition_matrix(max(quiet_hz/fs, 1e-9), p_burst, dt*fs/P, dt/tau)
    cuts = np.arange(phi % P, n, P); edges = np.r_[0, cuts, n]; cen = []
    d = np.diff(iq); s2 = (0.5*(np.median(np.abs(d.real))+np.median(np.abs(d.imag)))/(0.67449*np.sqrt(2)))**2
    for lo, hi in zip(edges[:-1], edges[1:]):
        a, b = lo+W, hi-W
        if b-a < 2*m+20: cen.append(None); continue
        fit, _ = smoothed_fit(iq[a:b], m); cen.append((fit.origin+fit.direction*fit.mu_a, fit.origin+fit.direction*fit.mu_b))
    total = np.zeros(2*W+1); used = 0
    for k, c in enumerate(cuts):
        pre, post = cen[k], cen[k+1]
        if pre is None or post is None or c-W-1 < 0 or c+W+1 > n: continue
        prof = switch_profile(iq[c-W-1:c+W+1], pre, post, s2, trans)    # switch at j = 1..2W+1 -> offsets -W..+W
        total += prof[:2*W+1]-prof[:2*W+1].max(); used += 1
    return (phi-W+int(np.argmax(total))) % P if used else phi, used


def estimate_phase_hmm2(iq, fs, period_s, K=3, W_s=6e-3):
    """'hmm' estimate: the K best coarse minima, each refined by :func:`hmm_fine` within +-W_s, the winner by
    :func:`record_loglik`. Returns (first onset in [0, period) [s], dict of candidates and their log-likelihoods)."""
    P = int(round(period_s*fs)); ys, off, s2 = prepare(iq); step = 100; grid = np.arange(0, P, step)
    cost = np.array([profile_cost(ys[::step], off[::step], g, P, 200) for g in grid])
    loc = [i for i in range(cost.size) if cost[i] <= cost[i-1] and cost[i] <= cost[(i+1) % cost.size]]
    loc.sort(key=lambda i: cost[i]); picked = []; sep = int(round(2*W_s*fs))
    for i in loc:
        g = int(grid[i])
        if all(min(abs(g-q), P-abs(g-q)) >= sep for q in picked): picked.append(g)
        if len(picked) == K: break
    cands = [hmm_fine(iq, fs, g, P, W_s)[0] for g in picked]
    L = [record_loglik(iq, fs, c, P) for c in cands] if len(cands) > 1 else [0.]
    return ((cands[int(np.argmax(L))]-0.5)/fs) % (P/fs), dict(candidates=cands, loglik=L)


def estimate_led_phase(iq, sample_rate_hz, period_s, method='hmm'):
    """Estimated first LED onset in [0, period) (seconds) and a diagnostic.

    method='hmm' (default): the 3 best coarse candidates, each refined by the exact HMM change-point profile
    within +-3 ms of every onset; the winner by whole-record HMM log-likelihood. method='v1': the coarse
    candidate refined by the memoryless mixture change-point likelihood.
    """
    iq = np.asarray(iq)
    if method == 'v1':
        e = estimate(iq, sample_rate_hz, period_s); return e['fp'], e['depth']
    if method == 'hmm':
        return estimate_phase_hmm2(iq, sample_rate_hz, period_s, K=3, W_s=3e-3)
    raise ValueError("method must be 'hmm' or 'v1'")
