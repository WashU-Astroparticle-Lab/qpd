"""Blind learning of the two-branch emission model from an I/Q trace.

Nothing here uses the device parameters. The model is recovered from the
trace's own geometry, which is fixed by three facts about dispersive parity
readout under an offset-charge ramp:

1. The odd-parity dispersive shift is the even-parity one displaced by half
   a Cooper pair, and chi_even is period-1 and even in n_g. So the *signed*
   branch splitting is periodic in n_g with period 1, passing through zero
   twice per period (at the parity-blind points n_g = 0.25 mod 0.5). Under a
   linear ramp this becomes a periodic function of time with fold period
   ``P = 0.5 / slope``.
2. Because the drive is parked far off resonance (detuning 12.94 MHz against
   a linewidth of 243 kHz, so |Delta| = 53 kappa), the parity-dependent part
   of S21 reduces to one fixed complex vector scaled by a *real* factor
   ~ 1/Delta. The two branch means therefore differ along an essentially
   **fixed line** in the I/Q plane -- the difference vector reverses at each
   blind point but does not rotate -- and the whole problem projects onto one
   real axis. Note the naive argument (that the splitting is small compared
   with the linewidth) is false here: the splitting is 22 kappa. Measured
   direction spread over n_g in [0, 0.5] is 0.42 mrad.
3. The common mode barely moves (a few percent of the noise), so the parity
   signal lives in the trace's *spread*, not its mean. That is what makes the
   period recoverable: the projected variance is ``sigma^2 + h(t)^2 / 4``,
   modulated by a factor of a few over one cycle.

The resulting model is ``mu_even/odd(t) = c(t) +/- h(t)/2`` on the projected
axis, with ``h`` signed and alternating between consecutive half-cycles. The
overall sign of ``h`` is not identifiable, but it corresponds to a global
even/odd relabelling and leaves flip *timing* untouched.

Two deterministic nuisances ride on top of the alternation, and both are
carried here rather than cleaned up after the fact:

* ``phase_offset`` -- a piecewise-constant shift of the ramp phase, which is
  how an offset-charge jump acts (see :mod:`.segment`).
* ``reset_period`` / ``reset_phase`` -- a sawtooth ramp whose span is an odd
  number of half Cooper pairs re-labels the two branches at every reset, so
  the sign schedule needs one extra flip per ramp period (see :mod:`.ramp`).
  At a 500 Hz ramp these resets outnumber real tunnelling events by more than
  an order of magnitude, so they have to live in the model.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = [
    "validate_trace",
    "EmissionModel",
    "estimate_direction",
    "estimate_fold_period",
    "learn_emission_model",
]


@dataclass
class EmissionModel:
    """Projected two-branch emission model learned from the trace.

    The folded profiles are kept alongside the evaluated per-sample arrays so
    that the ramp phase can be re-aligned, or a reset schedule installed,
    without relearning anything.
    """

    direction: complex  # unit vector of the discrimination axis in I/Q
    origin: complex  # I/Q point the projection is measured from
    sigma: float  # per-quadrature noise std
    common_mode: np.ndarray  # (n,) c(t) on the projected axis
    separation: np.ndarray  # (n,) signed h(t); |h| = branch splitting
    fold_period: float | None  # [s] period of |h|; None if no ramp detected
    common_profile: np.ndarray | None = None  # (n_bins,) c vs folded phase
    magnitude_profile: np.ndarray | None = None  # (n_bins,) |h| vs phase
    blind_phase: float = 0.0  # folded phase of the parity-blind point
    phase_offset: np.ndarray | None = None  # (n,) ramp-phase offset [cycles]
    reset_period: float | None = None  # [s] sawtooth ramp period
    reset_phase: float = 0.0  # [s] time of the first reset within a period
    diagnostics: dict = field(default_factory=dict)

    def project(self, iq: np.ndarray) -> np.ndarray:
        """Map complex I/Q onto the real discrimination axis."""
        return np.real(np.conj(self.direction) * (np.asarray(iq) - self.origin))

    def branch_means(self) -> tuple[np.ndarray, np.ndarray]:
        return (self.common_mode + 0.5 * self.separation,
                self.common_mode - 0.5 * self.separation)

    @property
    def contrast(self) -> np.ndarray:
        """Per-sample branch separation in units of the noise std."""
        return np.abs(self.separation) / self.sigma

    def reset_sign(self, t: np.ndarray) -> np.ndarray:
        """Extra branch-label sign from the ramp resets before each time.

        A sawtooth spanning an odd number of half Cooper pairs maps each
        branch onto the other every time it resets. That is observationally
        identical to a parity flip, so it cannot be told apart event by event
        -- but it is *exactly* periodic, and folding it into the sign schedule
        removes it from the decoder's output entirely.
        """
        if self.reset_period is None:
            return np.ones(np.size(t))
        k = np.floor((np.asarray(t, dtype=float) - self.reset_phase)
                     / self.reset_period)
        return np.where(np.mod(k, 2.0) == 0, 1.0, -1.0)

    def evaluate(
        self, t: np.ndarray, phase_offset: float | np.ndarray = 0.0
    ) -> tuple[np.ndarray, np.ndarray]:
        """Common mode and signed splitting at times ``t``.

        ``phase_offset`` shifts the ramp phase in units of the fold period,
        which is how an offset-charge jump acts on this model: a jump of
        ``delta`` in n_g moves the phase by ``2*delta`` (one fold period is
        half a Cooper pair).
        """
        t = np.asarray(t, dtype=float)
        if self.fold_period is None or self.magnitude_profile is None:
            return (np.full(t.size, float(np.mean(self.common_mode))),
                    np.full(t.size, float(np.mean(self.separation))))
        n_bins = self.magnitude_profile.size
        u = t / self.fold_period - phase_offset
        phase = np.mod(u, 1.0)
        idx = np.minimum((phase * n_bins).astype(np.intp), n_bins - 1)
        magnitude = self.magnitude_profile[idx]
        common = (self.common_profile[idx] if self.common_profile is not None
                  else np.zeros_like(magnitude))
        cycles = np.floor(u - self.blind_phase)
        sign = np.where(np.mod(cycles, 2.0) == 0, 1.0, -1.0)
        return common, sign * self.reset_sign(t) * magnitude

    def refresh(self, t: np.ndarray) -> "EmissionModel":
        """Re-evaluate the per-sample arrays from the profiles in place."""
        offset = 0.0 if self.phase_offset is None else self.phase_offset
        self.common_mode, self.separation = self.evaluate(t, offset)
        return self

    def with_reset(self, t: np.ndarray, period: float,
                   phase: float) -> "EmissionModel":
        """Copy carrying a ramp-reset schedule, with arrays refreshed."""
        out = EmissionModel(**{**self.__dict__,
                              "reset_period": float(period),
                              "reset_phase": float(phase)})
        out.diagnostics = {**self.diagnostics}
        return out.refresh(t)

    def with_phase_offset(self, t: np.ndarray,
                          offset: np.ndarray) -> "EmissionModel":
        """Copy carrying a per-sample ramp-phase offset, arrays refreshed."""
        out = EmissionModel(**{**self.__dict__,
                               "phase_offset": np.asarray(offset, dtype=float)})
        out.diagnostics = {**self.diagnostics}
        return out.refresh(t)


def validate_trace(iq: np.ndarray, sample_rate: float,
                   min_samples: int = 3) -> np.ndarray:
    """Reject traces that would otherwise decode into fabricated events.

    A single non-finite sample is enough to poison everything downstream: the
    projection origin is the trace mean, so one NaN makes every emission NaN,
    and in :func:`~qpd.reconstruction.hmm.viterbi` a NaN comparison is always
    False, so the backtrace alternates state at every step and returns one
    "flip" per sample. Nothing about that failure looks abnormal in the
    diagnostics, so it has to be caught at the door. Real DAQ traces do contain
    dropouts.
    """
    iq = np.asarray(iq)
    if iq.ndim != 1 or iq.size < min_samples:
        raise ValueError(
            f"iq must be a 1-D trace with at least {min_samples} samples; "
            f"got shape {iq.shape}")
    if not np.isfinite(sample_rate) or sample_rate <= 0:
        raise ValueError(f"sample_rate must be positive and finite; "
                         f"got {sample_rate!r}")
    if not np.all(np.isfinite(iq)):
        bad = int(np.count_nonzero(~np.isfinite(iq)))
        raise ValueError(
            f"iq contains {bad} non-finite sample(s) (NaN or inf). One is "
            "enough to make the decoder emit a flip at every sample; drop or "
            "interpolate them before reconstructing.")
    return iq


def estimate_direction(iq: np.ndarray) -> tuple[complex, complex, float]:
    """Principal axis of the I/Q cloud, its origin, and the noise sigma.

    The cloud is an isotropic Gaussian smeared along one line by the branch
    splitting, so the minor principal axis carries noise only -- its standard
    deviation *is* sigma. This is the cleanest blind noise estimate available
    here: it needs no assumption about the flip rate or the ramp.
    """
    iq = np.asarray(iq)
    origin = complex(np.mean(iq))
    d = iq - origin
    x, y = d.real, d.imag
    cov = np.array([[np.mean(x * x), np.mean(x * y)],
                    [np.mean(x * y), np.mean(y * y)]])
    vals, vecs = np.linalg.eigh(cov)  # ascending
    major = vecs[:, 1]
    sigma = float(np.sqrt(max(vals[0], 1e-30)))
    return complex(major[0], major[1]), origin, sigma


def _fold_contrast(x2: np.ndarray, t: np.ndarray, period: float,
                   n_bins: int, phase_offset: float | np.ndarray = 0.0) -> float:
    """Spread of the folded mean of x^2 -- maximal at the true fold period."""
    idx = np.mod(t / period - phase_offset, 1.0) * n_bins
    idx = np.minimum(idx.astype(np.intp), n_bins - 1)
    counts = np.bincount(idx, minlength=n_bins)
    sums = np.bincount(idx, weights=x2, minlength=n_bins)
    ok = counts > 0
    profile = sums[ok] / counts[ok]
    return float(np.var(profile))


def _maximise_fold_contrast(u: np.ndarray, t: np.ndarray, p0: float,
                            rel_span: float, n_bins: int, n_iter: int = 5,
                            phase_offset: float | np.ndarray = 0.0) -> float:
    """Zoom in on the fold period that maximises the folded contrast."""
    lo, hi = p0 * (1.0 - rel_span), p0 * (1.0 + rel_span)
    best = p0
    for _ in range(n_iter):
        trial = np.linspace(lo, hi, 41)
        objs = [_fold_contrast(u, t, p, n_bins, phase_offset) for p in trial]
        best = float(trial[int(np.argmax(objs))])
        span = (hi - lo) / 20.0
        lo, hi = best - span, best + span
    return best


def refine_fold_period(
    x: np.ndarray,
    dt: float,
    period: float,
    phase_offset: np.ndarray | float = 0.0,
    rel_span: float = 3e-5,
    n_bins: int = 64,
) -> float:
    """Polish the fold period on the whole trace at a known phase offset.

    Only usable once the offset-charge jumps have been located: with the phase
    steps removed the full trace folds coherently again, and the whole duration
    can be brought to bear on the period.
    """
    u = x * x
    u = u - u.mean()
    t = np.arange(x.size) * dt
    return _maximise_fold_contrast(u, t, period, rel_span, n_bins,
                                   phase_offset=phase_offset)


def estimate_fold_period(
    x: np.ndarray,
    dt: float,
    period_bounds: tuple[float, float] | None = None,
    n_bins: int = 64,
    significance: float = 20.0,
    min_cycles_per_window: float = 2000.0,
    max_windows: int = 8,
) -> tuple[float | None, dict]:
    """Recover the period of the branch-splitting envelope, blind.

    ``x`` is the projected trace. Its square has mean
    ``sigma^2 + h(t)^2 / 4``, so the ramp shows up as a strong single tone in
    the periodogram of ``x^2``. The periodogram peak is located first, then
    refined by maximising the folded contrast, which pins the period far more
    tightly than the ``1/T`` frequency resolution -- necessary because the
    phase error accumulates over every cycle of the trace, and a fast ramp
    packs tens of thousands of cycles into a few seconds.

    The refinement is then repeated on each of several equal windows and the
    **median** taken. An offset-charge jump puts a phase step in the middle of
    the trace, and a whole-trace fold objective will happily trade that step
    against a small period error: on the reference scenario a single jump biases
    the whole-trace estimate by 3e-5, which is three quarters of a cycle of
    accumulated phase by the end. A jump falls in only one window, so the median
    is untouched by it.

    Returns ``(period, diagnostics)``; ``period`` is None when no significant
    modulation exists (a static n_g), which the caller should treat as a
    constant-separation model.
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    duration = n * dt
    u = x * x
    u = u - u.mean()

    spec = np.abs(np.fft.rfft(u))
    freqs = np.fft.rfftfreq(n, dt)
    lo_f, hi_f = 5.0 / duration, 0.4 / dt
    if period_bounds is not None:
        lo_f = max(lo_f, 1.0 / period_bounds[1])
        hi_f = min(hi_f, 1.0 / period_bounds[0])
    band = (freqs >= lo_f) & (freqs <= hi_f)
    if not np.any(band):
        return None, {"reason": "empty search band"}
    # argmax over the in-band slice only: a whole-array argmax returns the DC
    # bin (index 0, out of band) when the in-band spectrum is flat or all-NaN,
    # and 1/freqs[0] is a division by zero.
    band_idx = np.flatnonzero(band)
    k = int(band_idx[np.argmax(spec[band_idx])])
    denom = float(np.median(spec[band]))
    prominence = float(spec[k] / denom) if denom > 0 else 0.0
    if not np.isfinite(prominence) or freqs[k] <= 0:
        return None, {"reason": "no usable spectral peak", "peak_hz": float(freqs[k])}
    diag = {"peak_hz": float(freqs[k]), "prominence": prominence}
    if prominence < significance:
        diag["reason"] = "no significant modulation"
        return None, diag

    # Refine: scan the folded contrast around the periodogram peak. The
    # periodogram locates the tone to ~1/T in frequency; the fold objective is
    # sensitive to coherent phase across the trace and does far better.
    t = np.arange(n) * dt
    p_peak = 1.0 / float(freqs[k])
    rel = min(1.5 * p_peak / duration, 0.2)  # one FFT bin, in relative terms
    p_full = _maximise_fold_contrast(u, t, p_peak, rel, n_bins)

    n_cycles = duration / p_full
    n_win = int(np.clip(n_cycles // max(min_cycles_per_window, 1.0),
                        1, max_windows))
    if n_win <= 1:
        diag.update({"refined_period_s": p_full, "n_windows": 1})
        return p_full, diag

    width = n // n_win
    per_window = []
    for w in range(n_win):
        seg = u[w * width:(w + 1) * width]
        per_window.append(_maximise_fold_contrast(
            seg, t[:seg.size], p_full, 1e-4, n_bins))
    best_p = float(np.median(per_window))
    diag.update({
        "refined_period_s": best_p,
        "whole_trace_period_s": p_full,
        "n_windows": n_win,
        "window_period_spread": float(np.std(per_window)),
    })
    return best_p, diag


def _folded_profiles(x: np.ndarray, phase: np.ndarray, n_bins: int):
    """Per-phase-bin mean and variance of the projected trace."""
    idx = np.minimum((phase * n_bins).astype(np.intp), n_bins - 1)
    counts = np.bincount(idx, minlength=n_bins).astype(float)
    s1 = np.bincount(idx, weights=x, minlength=n_bins)
    s2 = np.bincount(idx, weights=x * x, minlength=n_bins)
    safe = np.maximum(counts, 1.0)
    mean = s1 / safe
    var = s2 / safe - mean ** 2
    return mean, var, counts


def _smooth_periodic(v: np.ndarray, width: int) -> np.ndarray:
    if width <= 1:
        return v
    k = np.ones(width) / width
    pad = np.concatenate([v[-width:], v, v[:width]])
    return np.convolve(pad, k, mode="same")[width:width + v.size]


def learn_emission_model(
    iq: np.ndarray,
    dt: float,
    n_bins: int = 96,
    smooth_bins: int = 3,
    period_bounds: tuple[float, float] | None = None,
    phase_offset: np.ndarray | None = None,
    fold_period: float | None = None,
    profile_slice: tuple[int, int] | None = None,
) -> EmissionModel:
    """Learn ``mu_even/odd(t)`` from the trace alone.

    Pipeline: principal axis and sigma -> fold period from the spectrum of the
    squared projection -> per-phase common mode and splitting magnitude from
    folded moments -> sign of the splitting from the alternation between
    consecutive half-cycles.

    ``phase_offset`` is a per-sample ramp-phase shift (in fold-period units)
    applied *before* folding. Passing the offsets found by :mod:`.segment`
    un-blends a trace containing offset-charge jumps: without it, samples from
    either side of a jump land in the same phase bin and the learned splitting
    profile is the average of two misaligned curves. ``fold_period`` skips the
    period search when it is already known.

    ``profile_slice`` restricts the samples used to build the folded profiles
    (the model is still evaluated over the whole trace). This breaks the
    chicken-and-egg at the start of the pipeline: the offsets cannot be found
    without a profile, but a whole-trace profile is already corrupted by the
    jump. It is worst for a jump near ``delta = 0.25``, which blends the profile
    with its own half-period shift -- the result is nearly period-1/2, so an
    offset of 0 becomes indistinguishable from 0.5 and the sign schedule is
    unrecoverable. A profile taken from one jump-free window has neither problem.
    """
    iq = np.asarray(iq)
    n = iq.size
    t = np.arange(n) * dt

    direction, origin, sigma = estimate_direction(iq)
    x = np.real(np.conj(direction) * (iq - origin))

    if fold_period is not None:
        period, period_diag = float(fold_period), {"supplied": True}
    else:
        period, period_diag = estimate_fold_period(
            x, dt, period_bounds=period_bounds
        )

    if period is None:
        # Static n_g: one splitting for the whole trace. The variance of the
        # projection is sigma^2 + h^2/4 with h constant.
        var = float(np.var(x))
        h = 2.0 * np.sqrt(max(var - sigma * sigma, 0.0))
        return EmissionModel(
            direction=direction, origin=origin, sigma=sigma,
            common_mode=np.full(n, float(np.mean(x))),
            separation=np.full(n, h),
            fold_period=None,
            diagnostics={"period": period_diag, "static_separation": h,
                         "max_contrast": h / sigma},
        )

    offset = 0.0 if phase_offset is None else np.asarray(phase_offset, float)
    phase = np.mod(t / period - offset, 1.0)
    sel = (slice(None) if profile_slice is None
           else slice(int(profile_slice[0]), int(profile_slice[1])))
    # The sampling grid only ever visits P/dt distinct phases when the fold
    # period is commensurate with the sample period, and then most of a fixed
    # bin grid is empty. Left unchecked the empty bins get a splitting of zero
    # and `argmin` anchors the whole sign schedule on one of them, which turns
    # the decode into noise -- silently, since every other diagnostic still
    # looks normal. Shrink the grid until every bin is populated.
    n_bins = int(n_bins)
    while n_bins > 8:
        mean_b, var_b, counts = _folded_profiles(x[sel], phase[sel], n_bins)
        if counts.min() > 0:
            break
        n_bins //= 2
    else:
        mean_b, var_b, counts = _folded_profiles(x[sel], phase[sel], n_bins)
    mean_b = _smooth_periodic(mean_b, smooth_bins)
    var_b = _smooth_periodic(var_b, smooth_bins)
    # var = sigma^2 + h^2/4 for a balanced two-branch mixture.
    mag_b = 2.0 * np.sqrt(np.maximum(var_b - sigma * sigma, 0.0))

    # The blind point is the phase of minimum splitting; the sign of h
    # reverses there, so it anchors the half-cycle alternation.
    blind_bin = int(np.argmin(mag_b))
    blind_phase = (blind_bin + 0.5) / n_bins

    model = EmissionModel(
        direction=direction,
        origin=origin,
        sigma=sigma,
        common_mode=np.empty(n),
        separation=np.empty(n),
        fold_period=period,
        common_profile=mean_b,
        magnitude_profile=mag_b,
        blind_phase=blind_phase,
        phase_offset=(None if phase_offset is None
                      else np.asarray(phase_offset, float)),
        diagnostics={
            "period": period_diag,
            "blind_phase": blind_phase,
            "max_contrast": float(mag_b.max() / sigma),
            "median_contrast": float(np.median(mag_b) / sigma),
            "phase_bin_counts_min": float(counts.min()),
            "n_phase_bins": int(n_bins),
        },
    )
    return model.refresh(t)
