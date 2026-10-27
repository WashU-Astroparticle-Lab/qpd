"""Waveform-only four-state HMM decoder for LED-pulsed parity data.

Every LED pulse redraws the offset charge, so the record is cut at each LED
onset and every window is decoded on its own: the static-bias parity x
quiet/burst HMM with emissions fitted to that window alone.

* With every option left at its default, this is the plain per-window decoder.
* The recommended settings (HMM-SF) are ``pin_quiet_rate=True,
  smooth_fit_s=2e-3, min_smoothed_contrast=2.5``.
* :func:`predict_led_hmm` needs the LED phase (``first_pulse_time_s``).
* :func:`predict_led_hmm_blind` estimates the phase from the waveform when only
  the LED period is known.
"""

from __future__ import annotations

import numpy as np

from qpd.reconstruction.burst_hmm import _PARITY, burst_transition_matrix, forward_backward_regime, viterbi_regime
from qpd.reconstruction.events import extract_flips
from qpd.reconstruction.hmm import decode_with_rate, gaussian_log_emissions
from qpd.reconstruction.static_bias import fit_two_blobs, reconstruct_parity_flips_static


def window_detectability(iq, p_flip_init):
    """Contrast times sqrt(decoded dwell in samples) of a fresh two-blob fit.

    The dwell comes from a two-state hard-EM decode (two rate iterations) on
    the window's own fit. If the window's charge lands near the parity-blind
    point, the decoder toggles on noise, the decoded dwell collapses, and this
    statistic collapses with it.
    """
    fit = fit_two_blobs(iq)
    n = iq.size
    _, p, _ = decode_with_rate(fit.project(iq), np.full(n, fit.mu_a), np.full(n, fit.mu_b),
                               fit.sigma, p_flip_init, 2)
    return float(fit.contrast * np.sqrt(1.0 / max(p, 1e-12)))


def _pinned_flips(x, mu_a, mu_b, sigma, sample_rate_hz, t0, p_quiet, p_burst, burst_rate_hz, burst_tau):
    """Burst-aware decode of a projected segment, with the quiet flip probability fixed (no seed, no hard-EM refinement)."""
    dt, n = 1.0 / sample_rate_hz, x.size
    log_emit4 = gaussian_log_emissions(x, np.full(n, mu_a), np.full(n, mu_b), sigma)[[0, 1, 0, 1]]
    trans = burst_transition_matrix(p_quiet, p_burst, burst_rate_hz * dt, dt / burst_tau)
    post4, _ = forward_backward_regime(log_emit4, trans)
    path = _PARITY[viterbi_regime(log_emit4, trans)].astype(np.int8)
    times, _ = extract_flips(path, post4[1] + post4[3], dt, t0=t0)
    return times


def _decode_pinned(iq, sample_rate_hz, t0, p_quiet, p_burst, burst_rate_hz, burst_tau):
    """Pinned-rate segment decode, with emissions from a two-blob fit on the raw segment."""
    fit = fit_two_blobs(iq)
    return _pinned_flips(fit.project(iq), fit.mu_a, fit.mu_b, fit.sigma, sample_rate_hz, t0,
                         p_quiet, p_burst, burst_rate_hz, burst_tau)


def smoothed_fit(iq, m):
    """Two-blob fit on an m-sample moving average of a segment, and its contrast |mu_a - mu_b| / sigma there.

    Averaging over a dwell-short, noise-long window leaves the branch centres in place but shrinks the
    noise by sqrt(m), so the centres are found where the raw samples overlap too much to fit.
    """
    fit = fit_two_blobs(np.convolve(iq, np.ones(m)/m, mode="valid"))
    return fit, float(abs(fit.mu_a-fit.mu_b)/max(fit.sigma, 1e-12))


def predict_led_hmm(iq, sample_rate_hz, *, first_pulse_time_s=0.050,
                    on_time_s=50e-6, off_time_s=0.300,
                    quiescent_rate_hz=20.0, burst_tau_s=3.7e-3,
                    dead_detectability=None, pin_quiet_rate=False,
                    smooth_fit_s=None, min_smoothed_contrast=None):
    """Decode parity x quiet/burst independently between known charge resets.

    The record is cut at the LED onsets ``first_pulse_time_s + k * (on_time_s
    + off_time_s)``, and each segment's emissions are fitted from that
    segment's waveform. No SNR, mean tunneling count, true charge or event
    labels are inputs. Spontaneous charge jumps within a segment are not
    corrected. Returns ``(flip_times, info)``.

    Options (all off by default):

    * ``dead_detectability`` reports nothing from a segment whose
      :func:`window_detectability` is below it. Segments whose charge lands
      near the parity-blind point otherwise produce mostly false flips.
    * ``pin_quiet_rate=True`` fixes the quiet-regime flip probability at
      ``quiescent_rate_hz / sample_rate_hz``. Without it, the rate is
      re-estimated per segment by hard EM, which inflates it on marginal
      segments.
    * ``smooth_fit_s`` (requires ``pin_quiet_rate=True``) fits each segment's
      emissions with :func:`smoothed_fit` over that many seconds, instead of on
      the raw samples. The raw samples are still decoded, projected on the
      smoothed fit's axis, with noise taken from their first differences.
      Segments shorter than twice the smoothing window report nothing.
    * ``min_smoothed_contrast`` reports a segment only if the smoothed fit's
      contrast reaches it.

    Recommended (HMM-SF): ``pin_quiet_rate=True, smooth_fit_s=2e-3,
    min_smoothed_contrast=2.5``.
    """
    iq = np.asarray(iq)
    if (iq.ndim != 1 or iq.size < 8 or not np.isfinite(iq).all()
            or not np.isfinite(sample_rate_hz) or sample_rate_hz <= 0):
        raise ValueError("Expected a finite one-dimensional IQ trace and positive sample rate")
    settings = [first_pulse_time_s, on_time_s, off_time_s,
                quiescent_rate_hz, burst_tau_s]
    if (not np.isfinite(settings).all() or first_pulse_time_s < 0
            or on_time_s <= 0 or off_time_s < 0 or quiescent_rate_hz < 0
            or burst_tau_s <= 0):
        raise ValueError("Invalid public acquisition settings")
    if dead_detectability is not None and not (np.isfinite(dead_detectability) and dead_detectability > 0):
        raise ValueError("dead_detectability must be None or positive and finite")
    if smooth_fit_s is not None and not (pin_quiet_rate and np.isfinite(smooth_fit_s) and smooth_fit_s > 0):
        raise ValueError("smooth_fit_s must be positive and finite, and needs pin_quiet_rate=True")
    if min_smoothed_contrast is not None and (smooth_fit_s is None or not np.isfinite(min_smoothed_contrast)):
        raise ValueError("min_smoothed_contrast must be finite and needs smooth_fit_s")
    m = None if smooth_fit_s is None else max(int(round(smooth_fit_s * sample_rate_hz)), 1)
    period = on_time_s + off_time_s
    p_flip_init = max(quiescent_rate_hz / sample_rate_hz, 1e-9)
    sample_times = np.arange(iq.size, dtype=float) / sample_rate_hz
    rises = np.arange(first_pulse_time_s, iq.size / sample_rate_hz, period)
    edges = np.unique(np.r_[0, np.searchsorted(sample_times, rises), iq.size])
    predictions = []
    degenerate_segments = 0
    skipped_samples = 0
    dead_segments = 0
    dead_samples = 0
    for start, stop in zip(edges[:-1], edges[1:]):
        if stop - start < 8:
            skipped_samples += int(stop - start)
            continue
        if (dead_detectability is not None
                and window_detectability(iq[start:stop], p_flip_init) < dead_detectability):
            dead_segments += 1
            dead_samples += int(stop - start)
            continue
        if m is not None:
            seg = iq[start:stop]
            fit, contrast = smoothed_fit(seg, m) if seg.size >= 2*m else (None, -np.inf)
            if fit is None or (min_smoothed_contrast is not None and contrast < min_smoothed_contrast):
                dead_segments += 1
                dead_samples += int(stop - start)
                continue
            x = fit.project(seg)
            predictions.append(_pinned_flips(x, fit.mu_a, fit.mu_b, np.median(np.abs(np.diff(x)))/(0.67449*np.sqrt(2)),
                                             sample_rate_hz, start / sample_rate_hz, p_flip_init, 0.3, 1.0 / period, burst_tau_s))
            continue
        if pin_quiet_rate:
            predictions.append(_decode_pinned(iq[start:stop], sample_rate_hz, start / sample_rate_hz, p_flip_init,
                                              0.3, 1.0 / period, burst_tau_s))
            continue
        result = reconstruct_parity_flips_static(
            iq[start:stop], sample_rate_hz, t0=start / sample_rate_hz,
            p_flip_init=p_flip_init,
            n_rate_iterations=4, min_confidence=0.0, segment_blocks=1,
            decoder="viterbi", burst_aware=True, p_burst=0.3,
            burst_rate_hz=1.0 / period, burst_tau=burst_tau_s,
        )
        # The degeneracy flag is reported as a diagnostic only, never used as a gate.
        predictions.append(result.flip_times)
        degenerate_segments += int(result.degenerate)
    events = np.concatenate(predictions) if predictions else np.empty(0)
    info = {"segments": len(edges) - 1,
            "degenerate_segments": degenerate_segments,
            "skipped_short_segment_samples": skipped_samples}
    if dead_detectability is not None or m is not None:
        info.update(dead_detectability=None if dead_detectability is None else float(dead_detectability),
                    dead_segments=dead_segments, dead_fraction=dead_samples / iq.size)
    if m is not None:
        info.update(smooth_fit_s=float(smooth_fit_s),
                    min_smoothed_contrast=None if min_smoothed_contrast is None else float(min_smoothed_contrast))
    if pin_quiet_rate:
        # The legacy degeneracy flag needs a fitted rate, so pinned segments do not report it.
        info.update(pin_quiet_rate=True, degenerate_segments=None)
    return events, info


def predict_led_hmm_blind(iq, sample_rate_hz, *, on_time_s=50e-6, off_time_s=0.300, early_s=100e-6, phase_method='hmm', **kwargs):
    """:func:`predict_led_hmm` without a known ``first_pulse_time_s``.

    The LED phase is estimated from the waveform, with the period known, by
    :func:`qpd.mlebench.led_phase.estimate_led_phase`. The record is then cut
    ``early_s`` before each estimated onset and decoded. An early cut is cheap,
    because the post-onset veto still covers the onset; a late one is not.

    ``info['first_pulse_time_s']`` is the schedule actually used; give that,
    not an external value, to any onset veto. ``phase_method`` is ``'hmm'``
    (default: an exact HMM change-point profile at every onset) or ``'v1'``
    (two-cluster profile plus a per-sample change-point likelihood). When the
    phase is known, call ``predict_led_hmm(..., first_pulse_time_s=...)``
    instead.
    """
    from qpd.mlebench.led_phase import estimate_led_phase
    period = on_time_s + off_time_s
    fp, diag = estimate_led_phase(iq, sample_rate_hz, period, method=phase_method)
    first = (fp - early_s) % period
    events, info = predict_led_hmm(iq, sample_rate_hz, first_pulse_time_s=first, on_time_s=on_time_s, off_time_s=off_time_s, **kwargs)
    info.update(first_pulse_time_s=first, estimated_onset_s=fp, phase_method=phase_method, phase_diagnostic=diag)
    return events, info
