"""Native-sample, no-veto adaptation of PR50 HMM-SF (smoothed fit).

Source: https://github.com/WashU-Astroparticle-Lab/qpd/pull/50
Revision: bed0f74422f23d8ca56da4f6a77381ff90307791, mlebench/led_hmm.py
and led_network/inputs.py under src/qpd/.
Preserves 2 ms emission fitting, contrast gate 2.5, fixed 20 Hz quiet prior,
and four-state transitions with p_burst=.3. No whitening, sample rebuilding,
resampling, phase estimation, or onset veto. The public nominal tau controls
geometric regime lifetime; EMG mu and sigma DO NOT enter this decoder. It is
not an EMG hazard HMM. Never pass per-record latent parameters here.

Explicit numerical deviation: floor raw difference-noise sigma at 1e-12 and
handle exactly constant segments without fitting a singular Gaussian. Moving
averages retain upstream np.convolve arithmetic for reproducibility.
"""
from dataclasses import dataclass

import numpy as np

from qpd.reconstruction.burst_hmm import (
    _PARITY, burst_transition_matrix, forward_backward_regime, viterbi_regime,
)
from qpd.reconstruction.events import extract_flips
from qpd.reconstruction.hmm import gaussian_log_emissions
from qpd.reconstruction.static_bias import fit_two_blobs

CHANNEL_NAMES = ("scaled_smoothed_contrast", "contrast_gate", "posterior_swing",
                 "valid", "segment_start", "baseline_event")
UPSTREAM_REVISION = "bed0f74422f23d8ca56da4f6a77381ff90307791"


@dataclass
class HMMSFResult:
    flip_times: np.ndarray
    channels: np.ndarray  # float32 [6, n], names above; no calibrated probabilities
    diagnostics: dict


def predict_hmmsf(iq, sample_rate_hz, *, first_pulse_time_s, on_time_s,
                  off_time_s, nominal_emg):
    """Return native-grid SF features and events using only public inputs.

    ``nominal_emg`` is the shared public mapping mu_s/sigma_s/tau_s, never the
    simulated record's varying EMG. Only tau_s affects this algorithm.
    Posterior swing is retained below the contrast gate as in PR50. Validity
    means a segment was fitted, not that its inferred parity is correct.
    Event raster marks the right endpoint of the decoded transition interval;
    precise refined times remain in flip_times. No edge joins separate fits.
    """
    z = np.asarray(iq)
    if z.ndim != 1 or z.size < 8 or not np.isfinite(z).all():
        raise ValueError("Expected a finite 1-D IQ trace with >=8 samples")
    fs, first, on, off = map(float, (sample_rate_hz, first_pulse_time_s, on_time_s, off_time_s))
    if not np.isfinite([fs, first, on, off]).all() or fs <= 0 or first < 0 or on <= 0 or off < 0 or not np.isfinite(on+off):
        raise ValueError("Invalid public acquisition settings")
    nominal = {k: float(nominal_emg[k]) for k in ("mu_s", "sigma_s", "tau_s")}
    if not np.isfinite(list(nominal.values())).all() or nominal["mu_s"] < 0 or nominal["sigma_s"] <= 0 or nominal["tau_s"] <= 0:
        raise ValueError("Invalid public nominal EMG")
    z = z.astype(np.complex128, copy=False)
    n, period = z.size, on+off
    rises = np.arange(first, n/fs, period)
    edges = np.unique(np.r_[0, np.searchsorted(np.arange(n)/fs, rises), n])
    m = max(int(round(.002*fs)), 1)
    out = np.zeros((len(CHANNEL_NAMES), n), np.float32)
    trans = burst_transition_matrix(max(20./fs, 1e-9), .3, 1./(period*fs), 1./(fs*nominal["tau_s"]))
    events = []
    info = dict(segments=len(edges)-1, short_segments=0, dead_segments=0,
                constant_segments=0, sigma_floor_segments=0, fitted_segments=0,
                sample_rate_hz=fs, nominal_emg=nominal, emg_mu_sigma_used=False,
                burst_tau_s=nominal["tau_s"], quiet_rate_hz=20., p_burst=.3,
                smooth_fit_s=.002, min_smoothed_contrast=2.5, veto=False,
                resampled=False, whitened=False, upstream_revision=UPSTREAM_REVISION,
                channel_names=list(CHANNEL_NAMES))
    for start, stop in zip(edges[:-1], edges[1:]):
        out[4, start] = 1
        seg = z[start:stop]
        if seg.size < max(8, 2*m):
            info["short_segments"] += 1
            continue
        if np.all(seg == seg[0]):
            info["constant_segments"] += 1
            info["dead_segments"] += 1
            continue
        smooth = np.convolve(seg, np.ones(m)/m, mode="valid")
        # The generic fitter needs eight smoothed samples (always true in the
        # specified 5--500 kHz domain, guarded for other public acquisitions).
        if smooth.size < 8:
            info["short_segments"] += 1
            continue
        fit = fit_two_blobs(smooth)
        c = float(abs(fit.mu_a-fit.mu_b)/max(fit.sigma, 1e-12))
        x = fit.project(seg)
        raw_sigma = float(np.median(np.abs(np.diff(x)))/(.67449*np.sqrt(2)))
        sigma = max(raw_sigma, 1e-12)
        info["sigma_floor_segments"] += int(raw_sigma < 1e-12)
        emit = gaussian_log_emissions(x, np.full(x.size, fit.mu_a), np.full(x.size, fit.mu_b), sigma)[[0, 1, 0, 1]]
        post, _ = forward_backward_regime(emit, trans)
        p = post[1]+post[3]
        path = _PARITY[viterbi_regime(emit, trans)].astype(np.int8)
        if not np.isfinite(c) or not np.isfinite(p).all():
            raise FloatingPointError("Nonfinite HMM-SF fit/posterior")
        passed = c >= 2.5
        out[0, start:stop] = c/(c+2.5)
        out[1, start:stop] = passed
        out[2, start:stop] = np.abs(np.diff(p, prepend=p[0]))
        out[3, start:stop] = 1
        info["fitted_segments"] += 1
        info["dead_segments"] += int(not passed)
        if passed:
            times, _ = extract_flips(path, p, 1./fs, t0=start/fs)
            events.append(times)
            out[5, start+np.flatnonzero(np.diff(path) != 0)+1] = 1
    return HMMSFResult(np.concatenate(events) if events else np.empty(0), out, info)
