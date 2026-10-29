"""Model inputs computed from the waveform and the LED schedule alone.

* :func:`hmm_channels` -- three per-sample summaries of the LED-window HMM posterior (``qpd.mlebench.led_hmm`` with
  default options): posterior swing, burst-regime probability and parity confidence.
* :func:`hmm_sf_channels` -- the same idea for HMM-SF (2 ms smoothed-fit emissions, quiet rate pinned at 20 Hz, contrast
  gate 2.5): smoothed-fit contrast, gate and posterior swing.
* :func:`record_inputs` -- the record's burst intensity ``mu_hat`` (mean flips per LED pulse, a grid maximum-likelihood
  estimate over the gated LED windows) and the canonical rotation of the I/Q plane: the parity axis is rotated to
  124 deg and the 180-deg sign is fixed by how the blob midpoint moves with the separation across windows.
* :func:`normalize_waveform` -- per-channel median / MAD-of-differences normalisation.

All parity summaries are sign-invariant, because the parity label of an LED window is arbitrary.
"""
import math

import numpy as np
from scipy.stats import exponnorm

from qpd.mlebench import led_hmm as SF

try:
    from numba import njit
except ImportError:  # pragma: no cover - pure-Python fallback (slow, same arithmetic)
    def njit(*a, **k):
        return (lambda f: f) if not (len(a) == 1 and callable(a[0])) else a[0]

FS = 1e5
ON_TIME_S, OFF_TIME_S = 50e-6, .300
PERIOD = ON_TIME_S+OFF_TIME_S
QUIESCENT_RATE_HZ, BURST_TAU_S, P_BURST = 20., 3.7e-3, .3
QRATE = 20.0
SMOOTH_FIT_S, SMOOTH_M = 2e-3, 200
GATE = MIN_CONTRAST = 2.5
MU_GRID = np.r_[0.0, np.geomspace(0.05, 60.0, 72)]
POP = float((0+1+np.sqrt(10)+10)/4)   # mu_hat when no LED window passes the gate
MU_SCALE = math.log1p(10.0)            # network feature = log1p(mu)/MU_SCALE
TARGET_AXIS_DEG = 124.0
V4_REF_DEG = -35.41553634478914

# Burst-delay law (ex-Gaussian, 1.2 ms + 0.4 ms Gaussian with a 3.7 ms tail) smoothed over 50 us: the hazard kernel of mu_hat
_grid = np.arange(0, 0.05, 1e-6)
_K = 3.7e-3/0.4e-3
_pdf = exponnorm.pdf(_grid, _K, loc=1.2e-3, scale=0.4e-3)
_pdf /= exponnorm.sf(0, _K, loc=1.2e-3, scale=0.4e-3)
DENS = np.convolve(_pdf, np.ones(50)/50)[:len(_grid)]
DGRID = _grid


def segment_edges(n, sample_rate_hz, first_pulse_time_s):
    times = np.arange(n, dtype=float)/sample_rate_hz
    rises = np.arange(first_pulse_time_s, n/sample_rate_hz, PERIOD)
    return np.unique(np.r_[0, np.searchsorted(times, rises), n])


def hmm_channels(iq, sample_rate_hz, first_pulse_time_s):
    """float16 [3, N]: |P_t - P_{t-1}|, P[burst regime], |2 P_t - 1| of the per-window four-state HMM (P = P[branch B]).

    Zero at each window's first sample (swing) and in windows shorter than 8 samples.
    """
    from qpd.reconstruction.static_bias import reconstruct_parity_flips_static
    iq = np.asarray(iq)
    if iq.ndim != 1 or iq.size < 8 or not np.isfinite(iq).all():
        raise ValueError('Expected a finite one-dimensional IQ trace')
    out = np.zeros((3, iq.size), dtype=np.float32)
    edges = segment_edges(iq.size, sample_rate_hz, first_pulse_time_s)
    for start, stop in zip(edges[:-1], edges[1:]):
        if stop-start < 8:
            continue
        r = reconstruct_parity_flips_static(
            iq[start:stop], sample_rate_hz, t0=start/sample_rate_hz, p_flip_init=max(QUIESCENT_RATE_HZ/sample_rate_hz, 1e-9),
            n_rate_iterations=4, min_confidence=0.0, segment_blocks=1, decoder='viterbi', burst_aware=True, p_burst=P_BURST,
            burst_rate_hz=1.0/PERIOD, burst_tau=BURST_TAU_S)
        p = np.asarray(r.posterior, dtype=float); b = np.asarray(r.diagnostics['burst_posterior'], dtype=float)
        if p.shape != (stop-start,) or b.shape != (stop-start,) or not (np.isfinite(p).all() and np.isfinite(b).all()):
            raise ValueError('HMM posterior shape/finiteness contract changed')
        swing = np.abs(np.diff(p, prepend=p[0]))
        out[0, start:stop], out[1, start:stop], out[2, start:stop] = swing, b, np.abs(2*p-1)
    return np.clip(out, 0., 1.).astype(np.float16)


def segment_posterior(seg, sample_rate_hz, m):
    """HMM-SF on one window up to the posterior: smoothed fit, projection, pinned four-state forward-backward."""
    fit, contrast = SF.smoothed_fit(seg, m)
    x = fit.project(seg); dt = 1.0/sample_rate_hz
    sigma = np.median(np.abs(np.diff(x)))/(0.67449*np.sqrt(2))
    log_emit4 = SF.gaussian_log_emissions(x, np.full(x.size, fit.mu_a), np.full(x.size, fit.mu_b), sigma)[[0, 1, 0, 1]]
    trans = SF.burst_transition_matrix(max(QUIESCENT_RATE_HZ/sample_rate_hz, 1e-9), P_BURST, dt/PERIOD, dt/BURST_TAU_S)
    post4, _ = SF.forward_backward_regime(log_emit4, trans)
    return contrast, np.asarray(post4[1]+post4[3], dtype=float)


def hmm_sf_channels(iq, sample_rate_hz, first_pulse_time_s):
    """float16 [3, N]: c/(c+2.5) of the window's smoothed-fit contrast c, the gate (c >= 2.5), and |P_t - P_{t-1}|.

    Windows shorter than twice the 2 ms smoothing window get zeros.
    """
    iq = np.asarray(iq)
    if iq.ndim != 1 or iq.size < 8 or not np.isfinite(iq).all():
        raise ValueError('Expected a finite one-dimensional IQ trace')
    m = max(int(round(SMOOTH_FIT_S*sample_rate_hz)), 1); out = np.zeros((3, iq.size), dtype=np.float32)
    edges = segment_edges(iq.size, sample_rate_hz, first_pulse_time_s)
    for start, stop in zip(edges[:-1], edges[1:]):
        if stop-start < max(8, 2*m):
            continue
        contrast, p = segment_posterior(iq[start:stop], sample_rate_hz, m)
        if p.shape != (stop-start,) or not np.isfinite(p).all() or not np.isfinite(contrast):
            raise ValueError('HMM-SF posterior shape/finiteness contract changed')
        out[0, start:stop] = contrast/(contrast+MIN_CONTRAST)
        out[1, start:stop] = float(contrast >= MIN_CONTRAST)
        out[2, start:stop] = np.abs(np.diff(p, prepend=p[0]))
    return np.clip(out, 0., 1.).astype(np.float16)


@njit(cache=True)
def fwd2_ll(x, mu_a, mu_b, sig, pf):
    """Log-likelihood of a projected window under a two-state telegraph with per-step flip probabilities pf."""
    n = x.size; inv = 0.5/(sig*sig)
    a0 = .5; a1 = .5; ll = 0.
    for k in range(n):
        if k > 0:
            p = pf[k-1]; b0 = a0*(1-p)+a1*p; b1 = a1*(1-p)+a0*p
        else:
            b0 = a0; b1 = a1
        la = -inv*(x[k]-mu_a)**2; lb = -inv*(x[k]-mu_b)**2
        m = max(la, lb); x0 = b0*np.exp(la-m); x1 = b1*np.exp(lb-m); s = x0+x1
        ll += np.log(s)+m; a0 = x0/s; a1 = x1/s
    return ll


@njit(cache=True)
def ll_grid(x, mu_a, mu_b, sig, kern, mus, q):
    """fwd2_ll on a grid of burst intensities: hazard q + mu * kern(t - onset)."""
    n = x.size; out = np.empty(mus.size); pf = np.empty(max(n-1, 1)); dt = 1e-5
    nk = kern.size
    for g in range(mus.size):
        mu = mus[g]
        for k in range(n-1):
            h = q+(mu*kern[k] if k < nk else 0.)
            pf[k] = (1-np.exp(-2*h*dt))/2
        out[g] = fwd2_ll(x, mu_a, mu_b, sig, pf)
    return out


def segment_kernel(start, onset_s, n):
    k0 = int(np.ceil(onset_s*FS)); t = (np.arange(start, min(start+n, k0+5000)))/FS
    ker = np.zeros(min(n, max(0, k0+5000-start)))
    sel = np.arange(start, start+ker.size) >= k0
    ker[sel] = np.interp(t[sel]-onset_s, DGRID, DENS)
    return ker


def mle_from_ll(lls):
    """Grid argmax of the summed log-likelihoods, refined by a parabola in log(1 + mu)."""
    if not lls:
        return float('nan')
    lg = np.log1p(MU_GRID); L = np.sum(lls, axis=0); i = int(np.argmax(L))
    if 0 < i < len(L)-1:
        x0, x1, x2 = lg[i-1:i+2]; y0, y1, y2 = L[i-1:i+2]
        den = (x0-x1)*(x0-x2)*(x1-x2)
        A = (x2*(y1-y0)+x1*(y0-y2)+x0*(y2-y1))/den; B = (x2*x2*(y0-y1)+x1*x1*(y2-y0)+x0*x0*(y1-y2))/den
        xm = -B/(2*A) if A < 0 else x1
        return float(np.expm1(np.clip(xm, x0, x2)))
    return float(MU_GRID[i])


def record_inputs(iq, first_pulse_time_s):
    """Returns dict(mu_hat, mu_missing, n_gated, axis_deg, rotation_deg, sign_flipped, v4_margin, segments)."""
    iq = np.asarray(iq).astype(np.complex128); n = iq.size
    times = np.arange(n, dtype=float)/FS
    rises = np.arange(first_pulse_time_s, n/FS, PERIOD)
    idx = np.searchsorted(times, rises)
    edges = np.unique(np.r_[0, idx, n])
    lls, mids, seps, wts, axis_acc = [], [], [], [], 0j
    for start, stop in zip(edges[:-1], edges[1:]):
        L = stop-start
        if L < max(8, 2*SMOOTH_M):
            continue
        seg = iq[start:stop]; fit, c = SF.smoothed_fit(seg, SMOOTH_M)
        ca, cb = fit.origin+fit.direction*fit.mu_a, fit.origin+fit.direction*fit.mu_b; d = ca-cb
        if abs(d) > 0:
            axis_acc += c**2*(d/abs(d))**2
        mids.append((ca+cb)/2); seps.append(abs(d)); wts.append(L*min(c**2, 100.))
        k = np.flatnonzero(idx == start)
        if not k.size or c < GATE:
            continue
        onset = float(rises[k[0]]); x = fit.project(seg)
        sig = max(float(np.median(np.abs(np.diff(x)))/(0.67449*np.sqrt(2))), 1e-12)
        ll = ll_grid(x.astype(np.float64), float(fit.mu_a), float(fit.mu_b), sig, segment_kernel(start, onset, L), MU_GRID, QRATE)
        lls.append((ll-ll[0]).astype(np.float64))
    mu = mle_from_ll(lls); missing = not np.isfinite(mu)
    axis = float(np.degrees(np.angle(axis_acc))/2 % 180)
    phi = (TARGET_AXIS_DEG-axis+90.) % 180.-90.
    v4 = 0j
    if len(mids) >= 2:
        m, s, w = np.array(mids), np.array(seps), np.array(wts, float)
        v4 = complex(np.sum(w*(s-np.sum(w*s)/w.sum())*(m-np.sum(w*m)/w.sum())))
    margin = float(np.real(np.conj(np.exp(1j*np.radians(V4_REF_DEG)))*v4*np.exp(1j*np.radians(phi))))
    flipped = margin < 0
    rotation = (phi+(180. if flipped else 0.)+180.) % 360.-180.
    return dict(mu_hat=POP if missing else mu, mu_missing=bool(missing), n_gated=len(lls), axis_deg=axis, rotation_deg=float(rotation),
                sign_flipped=bool(flipped), v4_margin=margin, segments=len(mids))


def rotate(iq, rotation_deg):
    return np.asarray(iq).astype(np.complex128)*np.exp(1j*np.radians(rotation_deg))


def mu_feature(mu):
    return float(np.log1p(max(float(mu), 0.0))/MU_SCALE)


def normalize_waveform(waveform):
    """Per-channel (x - median) / sigma, sigma = MAD of first differences / (0.6745 sqrt 2); float32 [2, N]."""
    x = np.asarray(waveform, dtype=np.float32)
    if x.ndim != 2 or x.shape[0] != 2 or x.shape[1] < 2 or not np.isfinite(x).all():
        raise ValueError("Expected finite two-channel waveform")
    center = np.median(x, axis=-1, keepdims=True)
    noise = np.median(np.abs(np.diff(x, axis=-1)), axis=-1, keepdims=True) / (.67448975 * np.sqrt(2))
    return ((x - center) / noise.clip(1e-12)).astype(np.float32)
