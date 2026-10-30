"""Per-record health and out-of-distribution flags for the LED network (report-only; nothing here feeds the network).

Per-window quantities come from the same HMM-SF 2 ms smoothed two-blob fits that :func:`inputs.record_inputs` uses
(windows between LED onsets, at least 400 samples). Hard flags mean the network output should not be trusted without
inspection (``LEDReconstructionResult.degenerate``); soft flags are informational. Thresholds were set on simulated
in-distribution records (SNR 0.5-5, 0-10 tunnelling events per pulse, 1/2/5 s), where hard flags never fire.
"""
import math

import numpy as np

from . import inputs as MI
from .network import PERIOD_S

FS = MI.FS
VETO_S = 300e-6
QUIET_AFTER_S = 0.030                       # quiet region: > 30 ms after an onset (burst kernel ~ e^{-t/3.7 ms})
STANDARD_N = (100000, 200000, 500000)
WHITEN_AUTO_THRESHOLD = 0.05

THRESHOLDS = {
    'sample_rate_mismatch': dict(rationale='The network, its input channels and mu_hat assume 100 kHz sampling. Hard.'),
    'period_mismatch': dict(tol_s=1e-9, rationale='The slot geometry, the input channels and mu_hat assume a 0.30005 s LED period. Hard.'),
    'nonstandard_length': dict(rationale='Trained on 1, 2 and 5 s records; other lengths are tiled the same way but untested. Soft.'),
    'first_pulse_outside_period': dict(rationale='Training first pulses were in [0, 0.30005) s; a later one implies an unlogged earlier pulse. Soft.'),
    'axis_coherence': dict(min=0.8, rationale='Per-window parity axes disagree, so the canonical rotation is ill-defined; the network '
                                               'fails for axis errors >= 60 deg. In-distribution >= 0.83. Hard.'),
    'raw_axis': dict(lo=120.0, hi=131.0, rationale='Raw parity axis outside the simulator frame (124.7 +- 1 deg). Expected on real '
                                                     'data with an arbitrary LO phase; the canonical rotation handles it. Soft, informational.'),
    'v4_sign': dict(min_abs_cos=0.5, min_strength=0.25, rationale='The 180-deg sign of the canonical frame is weakly determined '
                                                                   '(a wrong sign is a mild loss). Soft.'),
    'coloured_noise': dict(max_abs_rho1=0.05, rationale='Lag-1 autocorrelation of the noise-only minor-axis residual of the model input. '
                                                        'Every noise scale in the pipeline is the MAD of first differences, which is biased '
                                                        'low for coloured noise, so likelihoods become over-confident (false-positive floods). '
                                                        'Fires when coloured noise is left unwhitened or whitening failed. In-distribution <= 0.007. Hard.'),
    'prewhitened_input': dict(rationale='I/Q were AR(1)-pre-whitened before all preprocessing: noise is white again but parity steps '
                                        'shrink by sqrt((1-rho)/(1+rho)), so the effective SNR is lower. Soft, informational.'),
    'snr': dict(lo=0.3, rationale='Estimated per-sample SNR below the trained 0.5-5 (the estimate reads ~5% low). Soft.'),
    'snr_high': dict(hi=6.0, rationale='Estimated per-sample SNR above the trained 0.5-5. The network degrades there (crowded bursts '
                                       'first, then everything: F1 0.77 at SNR 16) while HMM-SF becomes near-perfect, so use HMM-SF '
                                       '(fallback_hmmsf=True does it automatically). Hard.'),
    'mu': dict(max=30.0, rationale='mu_hat (flips per pulse) far above the trained 0-10, or missing (no gated window). Soft.'),
    'no_gated_segment': dict(rationale='No window reaches smoothed-fit contrast 2.5: mu_hat falls back to the population mean and the '
                                       'rotation comes from ungated fits. Hard.'),
    'agreement': dict(min=0.4, min_events=10, rationale='Matched fraction (200 us) between network and HMM-SF flips in gated windows. '
                                                         'HMM-SF is the weaker decoder, so agreement falls with SNR. Soft.'),
    'agreement_high_snr': dict(min=0.6, min_snr=2.0, rationale='At SNR >= 2 the two decoders should agree; low agreement means the '
                                                                'network is off-manifold (in-distribution >= 0.87). Hard.'),
    'parity': dict(min=0.4, min_segments=5, rationale='Parity of the decoded flip count per window vs the change of the HMM-SF decided '
                                                       'parity state across it (random would give 0.5). Soft, low power.'),
}


def segments(n, first_pulse_time_s):
    """(start, stop, onset_s or None) exactly as :func:`inputs.record_inputs` splits a record."""
    times = np.arange(n, dtype=float)/FS
    rises = np.arange(first_pulse_time_s, n/FS, PERIOD_S)
    idx = np.searchsorted(times, rises)
    edges = np.unique(np.r_[0, idx, n])
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        k = np.flatnonzero(idx == a)
        out.append((int(a), int(b), float(rises[k[0]]) if k.size else None))
    return out


def onset_times(first_pulse_time_s, duration_s, period_s=PERIOD_S):
    """Scheduled LED onsets in [0, duration_s)."""
    count = max(int(np.ceil((duration_s-first_pulse_time_s)/period_s)), 0)
    onsets = first_pulse_time_s + period_s*np.arange(count)
    return onsets[onsets < duration_s]


def veto_mask(times, onsets, window_s=VETO_S):
    """True for times outside [onset, onset + window_s) of the latest onset at or before them.

    Each LED pulse resets the offset charge, which steps the I/Q level like a parity flip; detections in the first
    300 us after an onset are dropped.
    """
    times = np.asarray(times, dtype=float)
    k = np.searchsorted(onsets, times, side='right')-1
    since = np.where(k >= 0, times-onsets[np.maximum(k, 0)], np.inf)
    return since >= window_s


def _lag1(y):
    y = y-y.mean(); d = float(np.dot(y, y))
    return float(np.dot(y[:-1], y[1:])/d) if d > 0 else float('nan')


def segment_table(iq, first_pulse_time_s, with_posterior=True):
    rows = []
    for a, b, onset in segments(iq.size, first_pulse_time_s):
        L = b-a
        if L < max(8, 2*MI.SMOOTH_M):
            continue
        seg = iq[a:b]
        fit, c = MI.SF.smoothed_fit(seg, MI.SMOOTH_M)
        ca, cb = fit.origin+fit.direction*fit.mu_a, fit.origin+fit.direction*fit.mu_b; d = ca-cb
        x = fit.project(seg)
        y = np.imag(np.conj(fit.direction)*(seg-fit.origin))          # minor axis: noise only
        sig = max(float(np.median(np.abs(np.diff(x)))/(0.67449*np.sqrt(2))), 1e-12)
        sig_perp_mad = float(np.median(np.abs(np.diff(y)))/(0.67449*np.sqrt(2)))
        sig_perp_std = float(np.std(y))
        r = dict(start=a, stop=b, onset_s=onset, contrast=c, gated=bool(c >= MI.GATE), sep=float(abs(d)), sigma=sig,
                 snr=float(abs(d)/sig), mid=(ca+cb)/2, d=d, weight=float(L*min(c**2, 100.)),
                 rho1=_lag1(y), sigma_ratio=float(sig_perp_mad/sig_perp_std) if sig_perp_std > 0 else float('nan'),
                 fit=fit, x=x)
        if with_posterior and r['gated']:
            _, p = MI.segment_posterior(seg, FS, MI.SMOOTH_M)
            r['posterior'] = p
        rows.append(r)
    return rows


def noise_rho(rows):
    """Length-weighted mean over LED windows of the lag-1 autocorrelation of the minor-axis projection of the raw samples.

    The two parity branches differ only along the major axis, so the minor axis carries noise only. Assumes the same
    noise colour in I and Q.
    """
    rho = np.array([r['rho1'] for r in rows if np.isfinite(r['rho1'])], float)
    w = np.array([r['stop']-r['start'] for r in rows if np.isfinite(r['rho1'])], float)
    return float(np.sum(w*rho)/w.sum()) if w.size else float('nan')


def whiten(I, Q, rho):
    """AR(1) pre-whitening: y_0 = x_0, y_k = (x_k - rho x_{k-1})/sqrt(1-rho^2). Returns float64 arrays."""
    out = []
    for x in (I, Q):
        x = np.asarray(x).astype(np.float32).astype(np.float64); y = x.copy()
        y[1:] = (x[1:]-rho*x[:-1])/np.sqrt(1-rho**2); out.append(y)
    return out[0], out[1]


def axis_and_sign(rows):
    acc = 0j; tot = 0.
    for r in rows:
        if abs(r['d']) > 0:
            acc += r['contrast']**2*(r['d']/abs(r['d']))**2; tot += r['contrast']**2
    axis = float(np.degrees(np.angle(acc))/2 % 180); coherence = float(abs(acc)/tot) if tot > 0 else float('nan')
    phi = (MI.TARGET_AXIS_DEG-axis+90.) % 180.-90.
    out = dict(axis_deg=axis, axis_coherence=coherence, segments=len(rows))
    if len(rows) >= 2:
        m = np.array([r['mid'] for r in rows]); s = np.array([r['sep'] for r in rows]); w = np.array([r['weight'] for r in rows])
        ds, dm = s-np.sum(w*s)/w.sum(), m-np.sum(w*m)/w.sum()
        v4 = complex(np.sum(w*ds*dm))
        margin = float(np.real(np.conj(np.exp(1j*np.radians(MI.V4_REF_DEG)))*v4*np.exp(1j*np.radians(phi))))
        den = math.sqrt(float(np.sum(w*ds**2))*float(np.sum(w*np.abs(dm)**2)))
        out.update(v4_margin=margin, v4_cos=margin/abs(v4) if abs(v4) > 0 else 0.0, v4_strength=abs(v4)/den if den > 0 else 0.0)
    else:
        out.update(v4_margin=0.0, v4_cos=0.0, v4_strength=0.0)
    return out


def quiet_rate_count(times, first_pulse_time_s, n):
    """Quiet flip rate from decoded flips > 30 ms after an onset (biased low by detection inefficiency)."""
    dur = n/FS; ons = onset_times(first_pulse_time_s, dur); prev = np.r_[first_pulse_time_s-PERIOD_S, ons]
    expo = 0.
    for k, o in enumerate(prev):
        lo = max(o+QUIET_AFTER_S, 0.); hi = min(prev[k+1] if k+1 < prev.size else dur, dur)
        expo += max(hi-lo, 0.)
    t = np.asarray(times, float); k = np.searchsorted(prev, t, side='right')-1
    quiet = (t-prev[np.maximum(k, 0)]) > QUIET_AFTER_S
    return dict(r_q_count_hz=float(quiet.sum()/expo) if expo > 0 else float('nan'), quiet_exposure_s=float(expo), quiet_flips=int(quiet.sum()))


def agreement(a, b, tol=200e-6):
    from qpd.mlebench.grade import match_times
    a, b = np.sort(np.asarray(a, float)), np.sort(np.asarray(b, float))
    if a.size+b.size == 0:
        return dict(matched=0, n_a=0, n_b=0, agreement=1.0)
    m = match_times(a, b, tol)
    return dict(matched=int(m.tp), n_a=int(a.size), n_b=int(b.size), agreement=float(2*m.tp/(a.size+b.size)))


def in_gated(times, rows):
    t = np.asarray(times, float); keep = np.zeros(t.size, bool)
    for r in rows:
        if r['gated']:
            keep |= (t >= r['start']/FS) & (t < r['stop']/FS)
    return t[keep]


def parity_consistency(times, rows, first_pulse_time_s, decided=0.95):
    """Fraction of decided gated windows where the parity of the decoded flip count matches the change of the HMM-SF
    decided parity state between the window start (after the veto) and its end."""
    t = np.asarray(times, float); n_dec = n_ok = 0
    for r in rows:
        p = r.get('posterior')
        if p is None:
            continue
        i0 = int(math.ceil(VETO_S*FS)) if r['onset_s'] is not None else 0
        if p.size <= i0+1:
            continue
        ps, pe = p[i0], p[-1]
        if min(abs(2*ps-1), abs(2*pe-1)) < 2*decided-1:
            continue
        lo, hi = (r['start']+i0)/FS, r['stop']/FS
        k = int(((t >= lo) & (t < hi)).sum())
        n_dec += 1; n_ok += int((k % 2 == 1) == ((ps > .5) != (pe > .5)))
    return dict(parity_segments_decided=n_dec, parity_consistent=n_ok,
                parity_consistency=float(n_ok/n_dec) if n_dec else float('nan'))


def flag(name, value, fired, hard, threshold, rationale):
    v = value if not (isinstance(value, float) and not math.isfinite(value)) else None
    return dict(name=name, value=v, fired=bool(fired), hard=bool(hard), threshold=threshold, rationale=rationale)


def evaluate(*, n, sample_rate_hz, first_pulse_time_s, period_s, mu_inputs, rows, net_times, hmm_times, whitening=None):
    """Return (flags, stats) for one record."""
    T = THRESHOLDS; ax = axis_and_sign(rows); F = []
    gated = [r for r in rows if r['gated']]
    snrs = [r['snr'] for r in gated] or [r['snr'] for r in rows]
    snr_hat = float(np.max(snrs)) if snrs else float('nan')
    rho1 = noise_rho(rows)
    sr = np.array([r['sigma_ratio'] for r in rows if np.isfinite(r['sigma_ratio'])])
    sigma_ratio = float(np.median(sr)) if sr.size else float('nan')
    ag_all = agreement(net_times, hmm_times); ag_g = agreement(in_gated(net_times, rows), in_gated(hmm_times, rows))
    pc = parity_consistency(net_times, rows, first_pulse_time_s)
    t = T['sample_rate_mismatch']; F.append(flag('sample_rate_mismatch', float(sample_rate_hz), abs(sample_rate_hz-FS) > 1e-6, True, FS, t['rationale']))
    t = T['period_mismatch']; F.append(flag('period_mismatch', float(period_s), abs(period_s-PERIOD_S) > t['tol_s'], True, PERIOD_S, t['rationale']))
    t = T['nonstandard_length']; F.append(flag('nonstandard_length', n/FS, n not in STANDARD_N, False, [1., 2., 5.], t['rationale']))
    t = T['first_pulse_outside_period']; F.append(flag('first_pulse_outside_period', float(first_pulse_time_s),
                                                   not (0 <= first_pulse_time_s < PERIOD_S), False, [0., PERIOD_S], t['rationale']))
    t = T['axis_coherence']; F.append(flag('axis_coherence_low', ax['axis_coherence'], not (ax['axis_coherence'] >= t['min']), True, t['min'], t['rationale']))
    t = T['raw_axis']; F.append(flag('raw_axis_off_training_frame', ax['axis_deg'], not (t['lo'] <= ax['axis_deg'] <= t['hi']), False, [t['lo'], t['hi']], t['rationale']))
    t = T['v4_sign']; F.append(flag('v4_sign_margin_small', ax['v4_cos'], not (abs(ax['v4_cos']) >= t['min_abs_cos'] and ax['v4_strength'] >= t['min_strength']), False,
                                    dict(min_abs_cos=t['min_abs_cos'], min_strength=t['min_strength']), t['rationale']))
    wh = whitening or {}
    t = T['coloured_noise']; F.append(flag('coloured_noise', rho1, not (abs(rho1) <= t['max_abs_rho1']), True, t['max_abs_rho1'], t['rationale']))
    t = T['prewhitened_input']; F.append(flag('prewhitened_input', wh.get('rho_used'), bool(wh.get('applied')), False, None, t['rationale']))
    t = T['snr']; F.append(flag('snr_below_trained_range', snr_hat, not (snr_hat >= t['lo']), False, t['lo'], t['rationale']))
    t = T['snr_high']; F.append(flag('snr_above_trained_range', snr_hat, snr_hat > t['hi'], True, t['hi'], t['rationale']))
    t = T['mu']; F.append(flag('mu_hat_high_or_missing', float(mu_inputs['mu_hat']), mu_inputs['mu_missing'] or mu_inputs['mu_hat'] > t['max'], False, t['max'], t['rationale']))
    t = T['no_gated_segment']; F.append(flag('no_gated_segment', len(gated), len(gated) == 0, True, 1, t['rationale']))
    t = T['agreement']; lowc = ag_g['n_a']+ag_g['n_b'] >= t['min_events']
    F.append(flag('network_hmmsf_disagree', ag_g['agreement'], lowc and ag_g['agreement'] < t['min'], False, t['min'], t['rationale']))
    t = T['agreement_high_snr']; hs = lowc and snr_hat >= t['min_snr']
    F.append(flag('network_hmmsf_disagree_high_snr', ag_g['agreement'], hs and ag_g['agreement'] < t['min'], True, dict(min=t['min'], min_snr=t['min_snr']), t['rationale']))
    t = T['parity']; v = pc['parity_consistency']
    F.append(flag('parity_count_inconsistent', v, pc['parity_segments_decided'] >= t['min_segments'] and v < t['min'], False, t['min'], t['rationale']))
    stats = dict(axis_deg=ax['axis_deg'], axis_coherence=ax['axis_coherence'], v4_margin=ax['v4_margin'], v4_cos=ax['v4_cos'], v4_strength=ax['v4_strength'],
                 segments=len(rows), gated_segments=len(gated), snr_hat=snr_hat, snr_hat_per_segment=[round(s, 4) for s in snrs],
                 noise_lag1_autocorr=rho1, sigma_madiff_over_std=sigma_ratio, whitening=wh,
                 agreement_all=ag_all, agreement_gated=ag_g, **pc)
    return F, stats
