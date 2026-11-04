"""One-call LED flip reconstruction with the trained network, shaped like the other ``qpd.reconstruction`` entry points."""
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import health as H
from . import inputs as MI
FS, ON_TIME_S, OFF_TIME_S, PERIOD_S = MI.FS, MI.ON_TIME_S, MI.OFF_TIME_S, MI.PERIOD
DEFAULT_WEIGHTS = Path(__file__).resolve().parent/'led_flip_network.pt'
_MODELS = {}


@dataclass
class LEDReconstructionResult:
    """Reconstructed flip times for one LED-pulsed record, and how they were obtained."""

    flip_times: np.ndarray            # [s] on the t0 axis, sorted, after the 300 us post-onset veto
    confidence: np.ndarray            # [0, 1] network presence probability per flip (NaN when the HMM-SF fallback answered)
    decoder: str = 'network'          # 'network' or 'hmm_sf_fallback'
    network_flip_times: np.ndarray = field(default_factory=lambda: np.empty(0))   # [s] the network's flips (always)
    hmmsf_flip_times: np.ndarray | None = None   # [s] the HMM-SF baseline's flips (when diagnostics are on)
    n_samples: int = 0
    sample_rate: float = FS
    t0: float = 0.0
    first_pulse_time_s: float = float('nan')   # [s] on the t0 axis: the given trigger time, or the estimate
    phase_estimated: bool = False             # True when first_pulse_time_s was estimated from the waveform
    diagnostics: dict = field(default_factory=dict)

    @property
    def rate_hz(self) -> float:
        """Quiet-period flip rate [Hz] counted from the reconstructed flips (> 30 ms after an onset; biased low)."""
        return float(self.diagnostics.get('quiet_rate', {}).get('r_q_count_hz', float('nan')))

    @property
    def degenerate(self) -> bool:
        """True when a hard health / out-of-distribution flag fired: do not use the flips without inspection.

        Requires ``diagnostics=True`` (the default); without it nothing is checked and this is False.
        """
        return bool(self.diagnostics.get('hard_flags_fired'))

    @property
    def flags(self) -> list:
        """Names of every flag that fired (hard first, then soft)."""
        return list(self.diagnostics.get('hard_flags_fired', []))+list(self.diagnostics.get('soft_flags_fired', []))

    @property
    def mu_hat(self) -> float:
        """Burst intensity fed to the network (mean tunnelling events per LED pulse)."""
        return float(self.diagnostics.get('mu', {}).get('mu_hat', float('nan')))

    @property
    def branch(self) -> np.ndarray:
        """(n,) int8 parity relative to the first sample, toggled at every reconstructed flip.

        The network outputs flip events, not a per-sample posterior, so this is derived from ``flip_times``; the
        absolute parity of the first sample is unknown.
        """
        idx = np.round((self.flip_times-self.t0)*self.sample_rate).astype(np.int64)
        return (np.searchsorted(idx, np.arange(self.n_samples), side='right') % 2).astype(np.int8)


def _model(weights, threads):
    try:
        import torch
    except ImportError as error:
        raise ImportError('The LED flip network needs torch: pip install "qpd[led]" (reconstruct_parity_flips_hmmsf does not)') from error
    if threads is not None:
        torch.set_num_threads(int(threads))
    key = str(weights)
    if key not in _MODELS:
        from .network import load_network
        _MODELS[key] = load_network(weights)
    return _MODELS[key]


def _veto(times, first_pulse_time_s, duration_s, extra=None, period_s=PERIOD_S):
    times = np.asarray(times, float)
    keep = H.veto_mask(times, H.onset_times(float(first_pulse_time_s), float(duration_s), period_s))
    return times[keep], (np.asarray(extra)[keep] if extra is not None else None)


def _hmm_sf(iq, first_pulse_time_s, sample_rate_hz=FS, on_time_s=ON_TIME_S, off_time_s=OFF_TIME_S):
    return MI.SF.predict_led_hmm(iq, sample_rate_hz, first_pulse_time_s=float(first_pulse_time_s), on_time_s=on_time_s, off_time_s=off_time_s,
                                 pin_quiet_rate=True, smooth_fit_s=2e-3, min_smoothed_contrast=2.5)


def preprocess(I, Q, first_pulse_time_s, period_s=PERIOD_S):
    """The network inputs: (x [2, N] float32, hmm channels, HMM-SF channels, record inputs, complex trace)."""
    I, Q = np.asarray(I).astype(np.float32), np.asarray(Q).astype(np.float32)
    if I.ndim != 1 or I.shape != Q.shape or I.size < 8 or not (np.isfinite(I).all() and np.isfinite(Q).all()):
        raise ValueError('I and Q must be finite 1-D arrays of equal length (>= 8)')
    fp = float(first_pulse_time_s)
    iq = I.astype(float)+1j*Q.astype(float)
    v1 = MI.hmm_channels(iq, FS, fp, period_s)
    sf = MI.hmm_sf_channels(iq, FS, fp, period_s)
    rec = MI.record_inputs(iq, fp, period_s)
    rot = MI.rotate(iq, rec['rotation_deg'])
    x = MI.normalize_waveform(np.stack([rot.real, rot.imag]))
    return x, v1, sf, rec, iq


def infer(model, x, v1, sf, mu_hat, first_pulse_time_s, on_time_s=ON_TIME_S, off_time_s=OFF_TIME_S):
    """Times [s] (before the veto) and presence of every slot with presence >= 0.5."""
    import torch
    from .network import configure, forward
    n = x.shape[-1]
    with torch.no_grad():
        configure(model, n, float(first_pulse_time_s), on_time_s=on_time_s, off_time_s=off_time_s)
        model.configure_hmm(torch.from_numpy(v1)[None])
        model.configure_sf(torch.from_numpy(sf)[None])
        model.configure_mu(MI.mu_feature(mu_hat))
        out = forward(model, torch.from_numpy(x)[None])[0].numpy()
    keep = out[:, 1] >= .5
    return out[keep, 0].astype(float)*(n/FS), out[keep, 1].astype(float)


def _split_iq(iq):
    a = np.asarray(iq) if not isinstance(iq, (tuple, list)) else np.stack([np.asarray(iq[0]), np.asarray(iq[1])])
    if np.iscomplexobj(a):
        if a.ndim != 1:
            raise ValueError('complex iq must be 1-D')
        return a.real, a.imag
    if a.ndim == 2 and a.shape[0] == 2:
        return a[0], a[1]
    raise ValueError('iq must be a 1-D complex trace, an (I, Q) pair or a (2, n) array')


def _estimate_phase(I, Q, sample_rate, period):
    from qpd.mlebench.led_phase import estimate_led_phase
    iq_est = np.asarray(I).astype(np.float32).astype(float)+1j*np.asarray(Q).astype(np.float32).astype(float)
    fp, est = estimate_led_phase(iq_est, sample_rate, period, method='hmm')
    fp = float(fp); ll = np.sort(np.asarray(est.get('loglik', []), float))[::-1]
    return fp, dict(estimated=True, method="qpd.mlebench.led_phase.estimate_led_phase(method='hmm')", first_pulse_time_s=fp,
                    loglik_margin=float(ll[0]-ll[1]) if ll.size > 1 else None)


def _finish_resampled(out, sample_rate, n, fallback_hmmsf):
    """Report a record resampled to 100 kHz on its own grid, and apply the fallback now that every flag is known."""
    diag = out.diagnostics
    out.n_samples, out.sample_rate = n, sample_rate
    diag.update(n_samples=n, sample_rate_hz=sample_rate)
    if 'flags' in diag:
        diag.update(hard_flags_fired=[f['name'] for f in diag['flags'] if f['fired'] and f['hard']],
                    soft_flags_fired=[f['name'] for f in diag['flags'] if f['fired'] and not f['hard']])
        if diag['hard_flags_fired'] and fallback_hmmsf:
            out.flip_times, out.confidence, out.decoder = out.hmmsf_flip_times, np.full(out.hmmsf_flip_times.size, np.nan), 'hmm_sf_fallback'
            diag.update(decoder='hmm_sf_fallback', n_flips=int(out.flip_times.size))
    return out


def _whitening_coefficient(rho_raw, whiten, whiten_rho, wh):
    """The AR(1) coefficient to whiten with, or None; records the choice in ``wh``."""
    wh.update(rho_hat_raw=rho_raw, auto_threshold=H.WHITEN_AUTO_THRESHOLD)
    if whiten_rho is not None:
        wh['source'] = 'whiten_rho (given)'
        return float(whiten_rho)
    if whiten == 'on' or (whiten == 'auto' and np.isfinite(rho_raw) and abs(rho_raw) > H.WHITEN_AUTO_THRESHOLD):
        wh['source'] = 'rho_hat'
        return rho_raw
    return None


def _rho_slow(I, Q, sample_rate, first_pulse_time_s, period, rebuild_seed):
    """Lag-1 noise autocorrelation of a record below 100 kHz, at its own rate."""
    ya, yb = MI.rebuild_100khz(I, Q, sample_rate, rebuild_seed)
    rows = H.segment_table(ya.astype(np.float32).astype(float)+1j*yb.astype(np.float32).astype(float), first_pulse_time_s, with_posterior=False, period_s=period)
    return H.noise_rho_slow(I, Q, rows, sample_rate)


def _reconstruct_fast(I, Q, sample_rate, first_pulse_time_s, t0, on_time_s, off_time_s, diagnostics, whiten, whiten_rho, fallback_hmmsf,
                      weights, threads):
    """Above 100 kHz: average down to 100 kHz, which keeps the signal and lowers the noise, and report on the original grid."""
    n = int(np.asarray(I).size)
    yI, yQ = MI.average_100khz(I, Q, sample_rate)
    out = reconstruct_parity_flips_led(yI+1j*yQ, FS, first_pulse_time_s=first_pulse_time_s, t0=t0, on_time_s=on_time_s, off_time_s=off_time_s,
                                       diagnostics=diagnostics, whiten=whiten, whiten_rho=whiten_rho, weights=weights, threads=threads)
    out.diagnostics['resampled'] = dict(method='averaged to 100 kHz', n_samples=int(yI.size))
    return _finish_resampled(out, sample_rate, n, fallback_hmmsf)


def _reconstruct_slow(I, Q, sample_rate, first_pulse_time_s, t0, on_time_s, off_time_s, diagnostics, whiten, whiten_rho, fallback_hmmsf,
                      rebuild_seed, weights, threads):
    """Below 100 kHz: whiten at the record's rate, rebuild it to 100 kHz, run the network, and report on the record's grid.

    Noise colour is measured and whitened here because the rebuilt record hides it. The HMM-SF baseline and the
    fallback run on the record's own samples, not on the rebuilt noise.
    """
    t_start = time.perf_counter()
    shift = float(t0); period = on_time_s+off_time_s; n = int(np.asarray(I).size); dur = n/sample_rate
    estimated = first_pulse_time_s is None
    if estimated:
        # On the rebuilt record: the estimator's window lengths are set in samples for 100 kHz.
        fp, phase = _estimate_phase(*MI.rebuild_100khz(I, Q, sample_rate, rebuild_seed), FS, period)
    else:
        fp = float(first_pulse_time_s)-shift
        phase = dict(estimated=False, first_pulse_time_s=fp)

    wh = dict(mode=whiten, applied=False, measured_at_hz=sample_rate)
    rho_after = float('nan')
    if whiten != 'off' or whiten_rho is not None or diagnostics:
        rho_after = rho_raw = _rho_slow(I, Q, sample_rate, fp, period, rebuild_seed)
        rho_use = _whitening_coefficient(rho_raw, whiten, whiten_rho, wh)
        if rho_use is not None and np.isfinite(rho_use) and abs(rho_use) < 1:
            I, Q = H.whiten(I, Q, rho_use); wh.update(applied=True, rho_used=rho_use)
            rho_after = wh['rho_hat_after'] = _rho_slow(I, Q, sample_rate, fp, period, rebuild_seed)
    yI, yQ = MI.rebuild_100khz(I, Q, sample_rate, rebuild_seed)
    out = reconstruct_parity_flips_led(yI+1j*yQ, FS, first_pulse_time_s=fp+shift, t0=shift, on_time_s=on_time_s, off_time_s=off_time_s, diagnostics=diagnostics,
                                       whiten='off', weights=weights, threads=threads)
    diag = out.diagnostics
    diag.update(phase=phase, whitening=wh,
                options=dict(diagnostics=diagnostics, whiten=whiten, whiten_rho=whiten_rho, fallback_hmmsf=fallback_hmmsf, rebuild_seed=rebuild_seed),
                resampled=dict(method='rebuilt to 100 kHz', seed=rebuild_seed, n_samples=int(yI.size), snr_units='per 100 kHz sample (per-sample SNR / sqrt(100 kHz / sample_rate))'))
    out.phase_estimated = estimated
    if diagnostics:
        hmm_raw, hmm_info = _hmm_sf(np.asarray(I, float)+1j*np.asarray(Q, float), fp, sample_rate, on_time_s, off_time_s)
        hmm, _ = _veto(hmm_raw, fp, dur, period_s=period); hmm = np.sort(hmm)
        out.hmmsf_flip_times = hmm+shift
        diag['hmmsf'] = dict(n_flips=int(hmm.size), info=hmm_info)
        diag['stats'].update(noise_lag1_autocorr=rho_after, whitening=wh)
        for f in diag['flags']:
            if f['name'] == 'coloured_noise':
                t = H.THRESHOLDS['coloured_noise']
                f.update(value=rho_after, fired=not (abs(rho_after) <= t['max_abs_rho1']))
            elif f['name'] == 'prewhitened_input':
                f.update(value=wh.get('rho_used'), fired=wh['applied'])
    diag['seconds'] = time.perf_counter()-t_start
    return _finish_resampled(out, sample_rate, n, fallback_hmmsf)


def reconstruct_parity_flips_led(iq, sample_rate=FS, first_pulse_time_s=None, t0=0.0, *, on_time_s=ON_TIME_S, off_time_s=OFF_TIME_S,
                                 diagnostics=True, whiten='auto', whiten_rho=None, fallback_hmmsf=False,
                                 rebuild_seed=0, weights=None, threads=None) -> LEDReconstructionResult:
    """Recover parity-flip times from one **LED-pulsed** I/Q record with the trained network.

    The network was trained on LED pulses of 50 us on, 0.300 s off; give another schedule with ``on_time_s`` and
    ``off_time_s``. Each pulse redraws the offset charge and releases a burst of quasiparticles; the network sees the
    whole record and reports every flip, quiet-period and burst alike. The trigger time (``first_pulse_time_s``) is best
    given; without it the phase is estimated.

    Parameters
    ----------
    iq : complex array, or (I, Q) pair / (2, n) array
        One record, ``I + 1j*Q`` on a uniform grid. Trained on 1, 2 and 5 s records.
    sample_rate : float
        Sampling rate [Hz]; each sample should average its full sample period. The network itself runs at 100 kHz. A
        faster record is averaged down to it; a slower one is rebuilt up to it (:func:`inputs.rebuild_100khz`), and flips
        closer than one sample cannot be separated.
    first_pulse_time_s : float, optional
        Time of the first LED pulse onset [s], on the same axis as ``t0`` (the LED trigger). If omitted, it is
        estimated from the waveform with :func:`qpd.mlebench.led_phase.estimate_led_phase` (``method='hmm'``, the
        period known) and every input -- LED cues, slot layout, HMM channels and the veto -- uses the estimate.
        This costs accuracy, mostly at low SNR where the estimate is occasionally off by milliseconds; pass the
        trigger time whenever it is logged. The value used is in ``result.first_pulse_time_s``.
    t0 : float
        Time of the first sample [s]; flip times are reported on this axis.
    on_time_s, off_time_s : float
        LED on and off times [s]. Every input, the slot layout, the veto and HMM-SF follow them, with the same weights.
        Off the trained schedule the soft flag ``untrained_schedule`` fires. An on time above 100 us or a period outside
        10 ms to 1 s raises the hard flag ``untested_schedule``: from 1 ms on, HMM-SF is the better decoder.
    diagnostics : bool
        Run the health / out-of-distribution flags and the HMM-SF baseline (default). False is fastest, and then
        ``degenerate`` cannot fire.
    whiten : 'auto' | 'on' | 'off'
        AR(1) pre-whitening of I and Q before all preprocessing. 'auto' whitens only when the estimated lag-1
        noise autocorrelation exceeds 0.05 in magnitude, so white-noise records are processed exactly as in
        training. ``whiten_rho`` forces a known coefficient.
    fallback_hmmsf : bool
        When a hard flag fires, return the HMM-SF baseline's flips instead (``decoder == 'hmm_sf_fallback'``).
    rebuild_seed : int
        Seed of the block noise when a record below 100 kHz is rebuilt, so the output is reproducible.
    weights : path, optional
        Network weights (default: the bundled ``led_flip_network.pt``).
    threads : int, optional
        torch CPU threads (outputs do not depend on it).

    Returns
    -------
    LEDReconstructionResult
        ``flip_times`` and per-flip ``confidence``; ``degenerate`` / ``flags`` (health); ``rate_hz``; ``mu_hat``;
        ``first_pulse_time_s`` / ``phase_estimated``; ``branch`` (derived); ``hmmsf_flip_times``; and the full
        ``diagnostics`` report.

    Examples
    --------
    >>> rec = reconstruct_parity_flips_led(iq, 1e5, first_pulse_time_s=0.1234)   # trigger time logged
    >>> rec = reconstruct_parity_flips_led(iq, 1e5)                              # phase estimated from the waveform
    >>> rec = reconstruct_parity_flips_led(iq, 1e4, first_pulse_time_s=0.0231, on_time_s=1e-4, off_time_s=60e-3)
    >>> if not rec.degenerate:
    ...     plot_trace_with_flips(iq, 1e5, rec.flip_times, confidence=rec.confidence, window=(0.12, 0.14))
    """
    if whiten not in ('auto', 'on', 'off'):
        raise ValueError("whiten must be 'auto', 'on' or 'off'")
    if fallback_hmmsf and not diagnostics:
        raise ValueError('fallback_hmmsf needs diagnostics (the hard flags)')
    t_start = time.perf_counter()
    on_time_s, off_time_s = float(on_time_s), float(off_time_s)
    if not (on_time_s > 0 and off_time_s > 0):
        raise ValueError('on_time_s and off_time_s must be positive')
    I, Q = _split_iq(iq)
    fs = float(sample_rate)
    if not fs > 0:
        raise ValueError('sample_rate must be positive')
    if fs > FS+1e-6:
        return _reconstruct_fast(I, Q, fs, first_pulse_time_s, t0, on_time_s, off_time_s, diagnostics, whiten, whiten_rho, fallback_hmmsf, weights, threads)
    if fs < FS-1e-6:
        return _reconstruct_slow(I, Q, fs, first_pulse_time_s, t0, on_time_s, off_time_s, diagnostics, whiten, whiten_rho, fallback_hmmsf,
                                 rebuild_seed, weights, threads)
    shift = float(t0); period = on_time_s+off_time_s
    n = int(np.asarray(I).size); dur = n/fs
    estimated = first_pulse_time_s is None
    if estimated:
        fp, phase = _estimate_phase(I, Q, fs, period)
    else:
        fp = float(first_pulse_time_s)-shift
        phase = dict(estimated=False, first_pulse_time_s=fp)
    diag = dict(n_samples=n, duration_s=dur, sample_rate_hz=fs, first_pulse_time_s=fp, t0=shift, on_time_s=on_time_s, off_time_s=off_time_s, led_period_s=period, phase=phase,
                options=dict(diagnostics=diagnostics, whiten=whiten, whiten_rho=whiten_rho, fallback_hmmsf=fallback_hmmsf),
                weights=str(weights or DEFAULT_WEIGHTS), veto_s=H.VETO_S)

    def result(flips, conf, decoder, net, hmm):
        return LEDReconstructionResult(
            flip_times=np.asarray(flips, float)+shift, confidence=np.asarray(conf, float), decoder=decoder,
            network_flip_times=np.asarray(net, float)+shift, hmmsf_flip_times=None if hmm is None else np.asarray(hmm, float)+shift,
            n_samples=n, sample_rate=fs, t0=shift, first_pulse_time_s=fp+shift, phase_estimated=estimated, diagnostics=diag)

    # noise colour: estimate rho from the per-window fits; optionally pre-whiten I/Q
    rows = None; wh = dict(mode=whiten, applied=False)
    if whiten != 'off' or whiten_rho is not None or diagnostics:
        iq0 = np.asarray(I).astype(np.float32).astype(float)+1j*np.asarray(Q).astype(np.float32).astype(float)
        rows = H.segment_table(iq0, fp, with_posterior=diagnostics, period_s=period)
        rho_use = _whitening_coefficient(H.noise_rho(rows), whiten, whiten_rho, wh)
        if rho_use is not None and np.isfinite(rho_use) and abs(rho_use) < 1:
            I, Q = H.whiten(I, Q, rho_use); wh.update(applied=True, rho_used=rho_use)
            iq0 = np.asarray(I).astype(np.float32).astype(float)+1j*np.asarray(Q).astype(np.float32).astype(float)
            rows = H.segment_table(iq0, fp, with_posterior=diagnostics, period_s=period)
            wh['rho_hat_after'] = H.noise_rho(rows)
    diag['whitening'] = wh
    model = _model(weights or DEFAULT_WEIGHTS, threads)
    x, v1, sf, rec, iq_model = preprocess(I, Q, fp, period)
    before, presence = infer(model, x, v1, sf, rec['mu_hat'], fp, on_time_s, off_time_s)
    net, net_p = _veto(before, fp, dur, presence, period)
    order = np.argsort(net, kind='stable'); net, net_p = net[order], net_p[order].astype(float)
    diag.update(mu=dict(mu_hat=rec['mu_hat'], mu_missing=rec['mu_missing'], n_gated=rec['n_gated'], feature=MI.mu_feature(rec['mu_hat'])),
                canonical=dict(axis_deg=rec['axis_deg'], rotation_deg=rec['rotation_deg'], sign_flipped=rec['sign_flipped'],
                               v4_margin=rec['v4_margin'], segments=rec['segments']),
                n_flips_before_veto=int(before.size))
    decoder, flips, conf, hmm = 'network', net, net_p, None
    if diagnostics:
        hmm_raw, hmm_info = _hmm_sf(iq_model, fp, FS, on_time_s, off_time_s)
        hmm, _ = _veto(hmm_raw, fp, dur, period_s=period); hmm = np.sort(hmm)
        flags, stats = H.evaluate(n=n, sample_rate_hz=fs, first_pulse_time_s=fp, period_s=period, mu_inputs=rec, rows=rows,
                                  net_times=net, hmm_times=hmm, whitening=wh, on_time_s=on_time_s)
        hard = [f['name'] for f in flags if f['fired'] and f['hard']]
        soft = [f['name'] for f in flags if f['fired'] and not f['hard']]
        diag.update(stats=stats, flags=flags, hard_flags_fired=hard, soft_flags_fired=soft,
                    hmmsf=dict(n_flips=int(hmm.size), info=hmm_info), network_n_flips=int(net.size))
        if hard and fallback_hmmsf:
            decoder, flips, conf = 'hmm_sf_fallback', hmm, np.full(hmm.size, np.nan)
    diag['quiet_rate'] = H.quiet_rate_count(net, fp, n, period)
    diag.update(decoder=decoder, n_flips=int(flips.size), seconds=time.perf_counter()-t_start)
    return result(flips, conf, decoder, net, hmm)


def reconstruct_parity_flips_hmmsf(iq, sample_rate=FS, first_pulse_time_s=None, t0=0.0, *, on_time_s=ON_TIME_S, off_time_s=OFF_TIME_S,
                                   whiten='auto', whiten_rho=None, rebuild_seed=0) -> LEDReconstructionResult:
    """Recover parity-flip times from one LED-pulsed I/Q record with HMM-SF alone. Does not need torch.

    The same call as :func:`reconstruct_parity_flips_led`, and the same flips as its ``hmmsf_flip_times``: a per-window
    four-state HMM with 2 ms smoothed-fit emissions, the quiet rate pinned at 20 Hz and a contrast gate of 2.5
    (:func:`qpd.mlebench.led_hmm.predict_led_hmm`), then the 300 us post-onset veto. Windows that fail the gate report
    no flips, so at low SNR it misses more than the network does; on long LED on times and very short periods it holds
    up better.

    Parameters
    ----------
    iq, sample_rate, first_pulse_time_s, t0, on_time_s, off_time_s, whiten, whiten_rho, rebuild_seed
        As in :func:`reconstruct_parity_flips_led`. A record above 100 kHz is averaged down to 100 kHz. One below it is
        decoded on its own samples; ``rebuild_seed`` only affects the phase estimate and the noise-colour measurement.

    Returns
    -------
    LEDReconstructionResult
        ``decoder == 'hmm_sf'``, ``flip_times`` (also in ``hmmsf_flip_times``), ``confidence`` all NaN and no ``mu_hat``.
        ``degenerate`` covers only coloured noise left unwhitened and no window passing the gate; a record with few
        passing windows is not flagged, see ``diagnostics['hmmsf']['info']['dead_fraction']``.

    Examples
    --------
    >>> rec = reconstruct_parity_flips_hmmsf(iq, 1e4, first_pulse_time_s=0.0231, on_time_s=1e-4, off_time_s=60e-3)
    >>> rec = reconstruct_parity_flips_hmmsf(iq, 1e4, on_time_s=1e-4, off_time_s=60e-3)   # phase estimated from the waveform
    >>> rec.flip_times, rec.first_pulse_time_s, rec.phase_estimated
    """
    if whiten not in ('auto', 'on', 'off'):
        raise ValueError("whiten must be 'auto', 'on' or 'off'")
    on_time_s, off_time_s, fs = float(on_time_s), float(off_time_s), float(sample_rate)
    if not (on_time_s > 0 and off_time_s > 0):
        raise ValueError('on_time_s and off_time_s must be positive')
    if not fs > 0:
        raise ValueError('sample_rate must be positive')
    I, Q = _split_iq(iq)
    I, Q = np.asarray(I, float), np.asarray(Q, float)
    if I.size < 8 or not (np.isfinite(I).all() and np.isfinite(Q).all()):
        raise ValueError('I and Q must be finite 1-D arrays of equal length (>= 8)')
    n = I.size
    if fs > FS+1e-6:
        yI, yQ = MI.average_100khz(I, Q, fs)
        out = reconstruct_parity_flips_hmmsf(yI+1j*yQ, FS, first_pulse_time_s, t0, on_time_s=on_time_s, off_time_s=off_time_s, whiten=whiten, whiten_rho=whiten_rho)
        out.diagnostics['resampled'] = dict(method='averaged to 100 kHz', n_samples=int(yI.size))
        return _finish_resampled(out, fs, n, False)
    slow = fs < FS-1e-6
    t_start = time.perf_counter()
    shift = float(t0); period = on_time_s+off_time_s; dur = n/fs
    estimated = first_pulse_time_s is None
    if estimated:
        fp, phase = _estimate_phase(*MI.rebuild_100khz(I, Q, fs, rebuild_seed), FS, period) if slow else _estimate_phase(I, Q, fs, period)
    else:
        fp = float(first_pulse_time_s)-shift
        phase = dict(estimated=False, first_pulse_time_s=fp)

    def single(a, b):
        return a.astype(np.float32).astype(float)+1j*b.astype(np.float32).astype(float)

    def rho_of(a, b):
        return _rho_slow(a, b, fs, fp, period, rebuild_seed) if slow else H.noise_rho(H.segment_table(single(a, b), fp, with_posterior=False, period_s=period))

    wh = dict(mode=whiten, applied=False)
    rho_after = rho_raw = rho_of(I, Q)
    rho_use = _whitening_coefficient(rho_raw, whiten, whiten_rho, wh)
    if rho_use is not None and np.isfinite(rho_use) and abs(rho_use) < 1:
        I, Q = H.whiten(I, Q, rho_use); wh.update(applied=True, rho_used=rho_use)
        rho_after = wh['rho_hat_after'] = rho_of(I, Q)
    raw, info = _hmm_sf(I+1j*Q if slow else single(I, Q), fp, fs, on_time_s, off_time_s)
    flips, _ = _veto(raw, fp, dur, period_s=period); flips = np.sort(flips)
    T = H.THRESHOLDS
    flags = [H.flag('coloured_noise', rho_after, not (abs(rho_after) <= T['coloured_noise']['max_abs_rho1']), True, T['coloured_noise']['max_abs_rho1'], T['coloured_noise']['rationale']),
             H.flag('prewhitened_input', wh.get('rho_used'), wh['applied'], False, None, T['prewhitened_input']['rationale']),
             H.flag('no_gated_segment', info['segments']-info['dead_segments'], info['dead_segments'] >= info['segments'], True, 1,
                    'No window reaches smoothed-fit contrast 2.5, so HMM-SF reports no flips. Hard.')]
    diag = dict(n_samples=n, duration_s=dur, sample_rate_hz=fs, first_pulse_time_s=fp, t0=shift, on_time_s=on_time_s, off_time_s=off_time_s, led_period_s=period, phase=phase,
                options=dict(whiten=whiten, whiten_rho=whiten_rho, rebuild_seed=rebuild_seed), veto_s=H.VETO_S, whitening=wh, decoder='hmm_sf', n_flips=int(flips.size),
                hmmsf=dict(n_flips=int(flips.size), info=info), flags=flags, hard_flags_fired=[f['name'] for f in flags if f['fired'] and f['hard']],
                soft_flags_fired=[f['name'] for f in flags if f['fired'] and not f['hard']], quiet_rate=H.quiet_rate_count(flips, fp, int(round(dur*FS)), period),
                seconds=time.perf_counter()-t_start)
    return LEDReconstructionResult(flip_times=flips+shift, confidence=np.full(flips.size, np.nan), decoder='hmm_sf', hmmsf_flip_times=flips+shift, n_samples=n, sample_rate=fs,
                                   t0=shift, first_pulse_time_s=fp+shift, phase_estimated=estimated, diagnostics=diag)
