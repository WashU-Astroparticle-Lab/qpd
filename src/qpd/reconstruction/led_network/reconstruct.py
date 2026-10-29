"""One-call LED flip reconstruction with the trained network, shaped like the other ``qpd.reconstruction`` entry points."""
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import health as H
from . import inputs as MI
from .network import PERIOD_S

FS = MI.FS
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
    import torch
    if threads is not None:
        torch.set_num_threads(int(threads))
    key = str(weights)
    if key not in _MODELS:
        from .network import load_network
        _MODELS[key] = load_network(weights)
    return _MODELS[key]


def _veto(times, first_pulse_time_s, duration_s, extra=None):
    times = np.asarray(times, float)
    keep = H.veto_mask(times, H.onset_times(float(first_pulse_time_s), float(duration_s)))
    return times[keep], (np.asarray(extra)[keep] if extra is not None else None)


def _hmm_sf(iq, first_pulse_time_s, sample_rate_hz=FS, period_s=PERIOD_S):
    kw = {} if sample_rate_hz == FS else dict(off_time_s=period_s-50e-6)
    return MI.SF.predict_led_hmm(iq, sample_rate_hz, first_pulse_time_s=float(first_pulse_time_s), pin_quiet_rate=True,
                                 smooth_fit_s=2e-3, min_smoothed_contrast=2.5, **kw)


def preprocess(I, Q, first_pulse_time_s):
    """The network inputs: (x [2, N] float32, hmm channels, HMM-SF channels, record inputs, complex trace)."""
    I, Q = np.asarray(I).astype(np.float32), np.asarray(Q).astype(np.float32)
    if I.ndim != 1 or I.shape != Q.shape or I.size < 8 or not (np.isfinite(I).all() and np.isfinite(Q).all()):
        raise ValueError('I and Q must be finite 1-D arrays of equal length (>= 8)')
    fp = float(first_pulse_time_s)
    iq = I.astype(float)+1j*Q.astype(float)
    v1 = MI.hmm_channels(iq, FS, fp)
    sf = MI.hmm_sf_channels(iq, FS, fp)
    rec = MI.record_inputs(iq, fp)
    rot = MI.rotate(iq, rec['rotation_deg'])
    x = MI.normalize_waveform(np.stack([rot.real, rot.imag]))
    return x, v1, sf, rec, iq


def infer(model, x, v1, sf, mu_hat, first_pulse_time_s):
    """Times [s] (before the veto) and presence of every slot with presence >= 0.5."""
    import torch
    from .network import configure, forward
    n = x.shape[-1]
    with torch.no_grad():
        configure(model, n, float(first_pulse_time_s))
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


def reconstruct_parity_flips_led(iq, sample_rate=FS, first_pulse_time_s=None, t0=0.0, *, led_period=PERIOD_S,
                                 diagnostics=True, whiten='auto', whiten_rho=None, fallback_hmmsf=False,
                                 weights=None, threads=None) -> LEDReconstructionResult:
    """Recover parity-flip times from one **LED-pulsed** I/Q record with the trained network.

    The LED pulses must follow the trained schedule (50 us on, 0.300 s off). Each pulse redraws the offset charge and
    releases a burst of quasiparticles; the network sees the whole record and reports every flip, quiet-period and
    burst alike. The trigger time (``first_pulse_time_s``) is best given; without it the phase is estimated.

    Parameters
    ----------
    iq : complex array, or (I, Q) pair / (2, n) array
        One record, ``I + 1j*Q`` on a uniform grid. Trained on 1, 2 and 5 s records.
    sample_rate : float
        Sampling rate [Hz]. The network is 100 kHz only; anything else raises unless ``fallback_hmmsf=True``.
    first_pulse_time_s : float, optional
        Time of the first LED pulse onset [s], on the same axis as ``t0`` (the LED trigger). If omitted, it is
        estimated from the waveform with :func:`qpd.mlebench.led_phase.estimate_led_phase` (``method='hmm'``, the
        period known) and every input -- LED cues, slot layout, HMM channels and the veto -- uses the estimate.
        This costs accuracy, mostly at low SNR where the estimate is occasionally off by milliseconds; pass the
        trigger time whenever it is logged. The value used is in ``result.first_pulse_time_s``.
    t0 : float
        Time of the first sample [s]; flip times are reported on this axis.
    led_period : float
        LED period [s]; anything other than 0.30005 s raises a hard flag.
    diagnostics : bool
        Run the health / out-of-distribution flags and the HMM-SF baseline (default). False is fastest, and then
        ``degenerate`` cannot fire.
    whiten : 'auto' | 'on' | 'off'
        AR(1) pre-whitening of I and Q before all preprocessing. 'auto' whitens only when the estimated lag-1
        noise autocorrelation exceeds 0.05 in magnitude, so white-noise records are processed exactly as in
        training. ``whiten_rho`` forces a known coefficient.
    fallback_hmmsf : bool
        When a hard flag fires, return the HMM-SF baseline's flips instead (``decoder == 'hmm_sf_fallback'``).
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
    >>> if not rec.degenerate:
    ...     plot_trace_with_flips(iq, 1e5, rec.flip_times, confidence=rec.confidence, window=(0.12, 0.14))
    """
    if whiten not in ('auto', 'on', 'off'):
        raise ValueError("whiten must be 'auto', 'on' or 'off'")
    if fallback_hmmsf and not diagnostics:
        raise ValueError('fallback_hmmsf needs diagnostics (the hard flags)')
    t_start = time.perf_counter()
    I, Q = _split_iq(iq)
    shift = float(t0); fs, period = float(sample_rate), float(led_period)
    n = int(np.asarray(I).size); dur = n/fs
    estimated = first_pulse_time_s is None
    if estimated:
        from qpd.mlebench.led_phase import estimate_led_phase
        iq_est = np.asarray(I).astype(np.float32).astype(float)+1j*np.asarray(Q).astype(np.float32).astype(float)
        fp, est = estimate_led_phase(iq_est, fs, period, method='hmm')
        fp = float(fp); ll = np.sort(np.asarray(est.get('loglik', []), float))[::-1]
        phase = dict(estimated=True, method="qpd.mlebench.led_phase.estimate_led_phase(method='hmm')", first_pulse_time_s=fp,
                     loglik_margin=float(ll[0]-ll[1]) if ll.size > 1 else None)
    else:
        fp = float(first_pulse_time_s)-shift
        phase = dict(estimated=False, first_pulse_time_s=fp)
    diag = dict(n_samples=n, duration_s=dur, sample_rate_hz=fs, first_pulse_time_s=fp, t0=shift, led_period_s=period, phase=phase,
                options=dict(diagnostics=diagnostics, whiten=whiten, whiten_rho=whiten_rho, fallback_hmmsf=fallback_hmmsf),
                weights=str(weights or DEFAULT_WEIGHTS), veto_s=H.VETO_S)

    def result(flips, conf, decoder, net, hmm):
        return LEDReconstructionResult(
            flip_times=np.asarray(flips, float)+shift, confidence=np.asarray(conf, float), decoder=decoder,
            network_flip_times=np.asarray(net, float)+shift, hmmsf_flip_times=None if hmm is None else np.asarray(hmm, float)+shift,
            n_samples=n, sample_rate=fs, t0=shift, first_pulse_time_s=fp+shift, phase_estimated=estimated, diagnostics=diag)

    if abs(fs-FS) > 1e-6:
        if not fallback_hmmsf:
            raise ValueError(f'sample_rate={fs}: the network is 100 kHz only (use fallback_hmmsf=True for HMM-SF at this rate)')
        times, _ = _hmm_sf(np.asarray(I, float)+1j*np.asarray(Q, float), fp, fs, period)
        ons = np.arange(fp, dur, period)
        flips = np.sort(times[H.veto_mask(times, ons)] if ons.size else times)
        diag.update(decoder='hmm_sf_fallback', n_flips=int(flips.size), hard_flags_fired=['sample_rate_mismatch'], soft_flags_fired=[],
                    flags=[H.flag('sample_rate_mismatch', fs, True, True, FS, H.THRESHOLDS['sample_rate_mismatch']['rationale'])])
        return result(flips, np.full(flips.size, np.nan), 'hmm_sf_fallback', np.empty(0), flips)
    # noise colour: estimate rho from the per-window fits; optionally pre-whiten I/Q
    rows = None; wh = dict(mode=whiten, applied=False)
    if whiten != 'off' or whiten_rho is not None or diagnostics:
        iq0 = np.asarray(I).astype(np.float32).astype(float)+1j*np.asarray(Q).astype(np.float32).astype(float)
        rows = H.segment_table(iq0, fp, with_posterior=diagnostics)
        rho_raw = H.noise_rho(rows); wh.update(rho_hat_raw=rho_raw, auto_threshold=H.WHITEN_AUTO_THRESHOLD)
        if whiten_rho is not None:
            rho_use = float(whiten_rho); wh['source'] = 'whiten_rho (given)'
        elif whiten == 'on' or (whiten == 'auto' and np.isfinite(rho_raw) and abs(rho_raw) > H.WHITEN_AUTO_THRESHOLD):
            rho_use = rho_raw; wh['source'] = 'rho_hat'
        else:
            rho_use = None
        if rho_use is not None and np.isfinite(rho_use) and abs(rho_use) < 1:
            I, Q = H.whiten(I, Q, rho_use); wh.update(applied=True, rho_used=rho_use)
            iq0 = np.asarray(I).astype(np.float32).astype(float)+1j*np.asarray(Q).astype(np.float32).astype(float)
            rows = H.segment_table(iq0, fp, with_posterior=diagnostics)
            wh['rho_hat_after'] = H.noise_rho(rows)
    diag['whitening'] = wh
    model = _model(weights or DEFAULT_WEIGHTS, threads)
    x, v1, sf, rec, iq_model = preprocess(I, Q, fp)
    before, presence = infer(model, x, v1, sf, rec['mu_hat'], fp)
    net, net_p = _veto(before, fp, dur, presence)
    order = np.argsort(net, kind='stable'); net, net_p = net[order], net_p[order].astype(float)
    diag.update(mu=dict(mu_hat=rec['mu_hat'], mu_missing=rec['mu_missing'], n_gated=rec['n_gated'], feature=MI.mu_feature(rec['mu_hat'])),
                canonical=dict(axis_deg=rec['axis_deg'], rotation_deg=rec['rotation_deg'], sign_flipped=rec['sign_flipped'],
                               v4_margin=rec['v4_margin'], segments=rec['segments']),
                n_flips_before_veto=int(before.size))
    decoder, flips, conf, hmm = 'network', net, net_p, None
    if diagnostics:
        hmm_raw, hmm_info = _hmm_sf(iq_model, fp)
        hmm, _ = _veto(hmm_raw, fp, dur); hmm = np.sort(hmm)
        flags, stats = H.evaluate(n=n, sample_rate_hz=fs, first_pulse_time_s=fp, period_s=period, mu_inputs=rec, rows=rows,
                                  net_times=net, hmm_times=hmm, whitening=wh)
        hard = [f['name'] for f in flags if f['fired'] and f['hard']]
        soft = [f['name'] for f in flags if f['fired'] and not f['hard']]
        diag.update(stats=stats, flags=flags, hard_flags_fired=hard, soft_flags_fired=soft,
                    hmmsf=dict(n_flips=int(hmm.size), info=hmm_info), network_n_flips=int(net.size))
        if hard and fallback_hmmsf:
            decoder, flips, conf = 'hmm_sf_fallback', hmm, np.full(hmm.size, np.nan)
    diag['quiet_rate'] = H.quiet_rate_count(net, fp, n)
    diag.update(decoder=decoder, n_flips=int(flips.size), seconds=time.perf_counter()-t_start)
    return result(flips, conf, decoder, net, hmm)
