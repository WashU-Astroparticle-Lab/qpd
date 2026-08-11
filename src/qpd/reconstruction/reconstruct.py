"""End-to-end blind reconstruction of parity-flip timing from an I/Q trace."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .emission import (EmissionModel, learn_emission_model,
                       refine_fold_period, validate_trace)
from .events import extract_flips
from .hmm import (decode_with_rate, decoded_path, forward_backward,
                  gaussian_log_emissions)
from .ramp import ResetComb, find_reset_comb
from .segment import segment_and_realign

__all__ = ["ReconstructionResult", "reconstruct_parity_flips_ramped"]


@dataclass
class ReconstructionResult:
    """Reconstructed parity-flip times and the model they came from."""

    flip_times: np.ndarray  # [s] estimated flip times
    confidence: np.ndarray  # [0, 1] posterior swing at each flip
    posterior: np.ndarray  # (n,) P[branch 1 | trace]
    branch: np.ndarray  # (n,) int8 Viterbi branch labels
    emission: EmissionModel | None = None
    p_flip: float = 0.0
    reset_comb: ResetComb | None = None
    charge_jump_times: np.ndarray = field(
        default_factory=lambda: np.empty(0, float))
    diagnostics: dict = field(default_factory=dict)

    @property
    def rate_hz(self) -> float:
        """Estimated per-state tunnelling rate."""
        return self.diagnostics.get("rate_hz", float("nan"))

    @property
    def sample_fidelity(self) -> float:
        """Mean single-sample parity assignment fidelity across the sweep.

        The swept case differs from the fixed-bias one in a way that matters
        here: the contrast is not a constant. It sweeps, passing through zero
        at every parity-blind crossing, so the per-sample assignment fidelity
        is a *function of ramp phase* -- unity where the branches are far apart
        and exactly 1/2 at the crossings. This property averages
        ``1 - erfc(C_k / (2*sqrt(2))) / 2`` over the trace.

        Because that average includes the blind crossings, it is a
        substantially pessimistic summary of the reconstruction, which does not
        classify samples independently (see :attr:`decoded_fidelity`). Use
        :meth:`fidelity_vs_phase` to see where in the sweep the information
        actually lives.
        """
        from scipy.special import erfc
        c = np.asarray(self.contrast, dtype=float)
        if c.size == 0 or not np.all(np.isfinite(c)):
            return float("nan")
        return float(np.mean(1.0 - 0.5 * erfc(c / (2.0 * np.sqrt(2.0)))))

    @property
    def decoded_fidelity(self) -> float:
        """Expected fraction of samples assigned to the correct branch.

        Read off the posterior as ``1 - mean(min(g, 1-g))``. Same caveat as in
        the fixed-bias case: this measures self-consistency with the fitted
        model, not correctness, so it stays high when the model itself is wrong.
        Read it together with :attr:`degenerate`.
        """
        g = np.asarray(self.posterior, dtype=float)
        if g.size == 0:
            return float("nan")
        return float(1.0 - np.mean(np.minimum(g, 1.0 - g)))

    def fidelity_vs_phase(self, n_bins: int = 32):
        """Single-sample assignment fidelity as a function of ramp phase.

        Returns ``(phase, fidelity)`` with ``phase`` the fold phase in [0, 1).
        This is the honest picture for a swept bias: the parity information is
        concentrated away from the blind crossings, and this shows how much of
        the sweep is actually useful. Returns ``(None, None)`` when no ramp was
        detected.
        """
        from scipy.special import erfc
        if self.emission is None or self.emission.fold_period is None:
            return None, None
        dt = float(self.diagnostics.get("dt", 0.0))
        if dt <= 0:
            return None, None
        c = np.asarray(self.contrast, dtype=float)
        f = 1.0 - 0.5 * erfc(c / (2.0 * np.sqrt(2.0)))
        phase = np.mod(np.arange(c.size) * dt / self.emission.fold_period, 1.0)
        idx = np.minimum((phase * n_bins).astype(np.intp), n_bins - 1)
        counts = np.bincount(idx, minlength=n_bins)
        sums = np.bincount(idx, weights=f, minlength=n_bins)
        ok = counts > 0
        return ((np.arange(n_bins)[ok] + 0.5) / n_bins, sums[ok] / counts[ok])

    @property
    def degenerate(self) -> bool:
        """True when the output should not be trusted as a list of events.

        The ramped pipeline has two ways to fail that leave every ordinary
        diagnostic looking healthy, so it needs an explicit self-check:

        * the reset comb locked onto a *multiple* of the true ramp period, so
          the sign schedule flips on only every j-th reset and the rest are
          emitted as flips;
        * the fold period is commensurate with the sample period, so the phase
          profile is built on a starved grid;
        * the trace is simply too noisy, in which case the decoder segments
          noise and the branch assignment is no better than a coin flip. This
          one produces *no* periodic structure, so it needs the contrast test
          rather than the comb test (measured: median contrast 0.67, decoded
          fidelity 0.96, true accuracy 0.500).

        The test for the first is that **no periodic comb may survive in the
        output**: real tunnelling is Poisson, so if the reported flip times
        still pile up at a fixed phase, the reset handling did not do its job.
        This is rate-invariant, unlike a threshold on the decoded rate.
        """
        return bool(self.diagnostics.get("degenerate", False))

    @property
    def contrast(self) -> np.ndarray:
        """Per-sample branch separation in units of the noise std."""
        return (self.emission.contrast if self.emission is not None
                else np.empty(0))


def _scan_reset_phase(x, emission, t, ramp_period, p_flip,
                      n_coarse=24, n_refine=2, scan_samples=200_000):
    """Locate the ramp-reset phase by maximising the HMM likelihood.

    Used when the ramp *period* is known from the hardware but its phase is
    not, which is the usual experimental situation: a signal generator's
    frequency is known far better than this analysis needs, while its trigger
    offset relative to the digitiser is not.

    Pinning the period turns comb-finding from a two-parameter search over a
    noisy first-pass event list into a one-parameter search -- and this version
    does not use the event list at all, it asks the trace directly which reset
    phase best explains it. That is what makes it keep working where the blind
    search has already failed: at a per-sample contrast of 0.67 the blind path
    loses the comb entirely (F1 0.02) while this recovers the phase exactly
    (F1 0.93).

    The scan runs on a prefix of the trace, since the phase is one fixed offset
    and a fraction of a second already holds hundreds of ramp cycles.
    """
    m = min(int(scan_samples), x.size)
    xs, ts = x[:m], t[:m]
    # The model may carry a per-sample ramp-phase offset (from charge-jump
    # segmentation); it has to be sliced to the prefix too, or the broadcast
    # against `ts` fails.
    if emission.phase_offset is not None and emission.phase_offset.size != m:
        emission = EmissionModel(**{**emission.__dict__,
                                    "phase_offset": emission.phase_offset[:m],
                                    "common_mode": emission.common_mode[:m],
                                    "separation": emission.separation[:m]})
    lo, hi = 0.0, ramp_period
    best_phi, best_ll = 0.0, -np.inf
    n_grid = int(n_coarse)
    for it in range(int(n_refine) + 1):
        for phi in np.linspace(lo, hi, n_grid, endpoint=(it > 0)):
            trial = emission.with_reset(ts, ramp_period, float(phi))
            a, b = trial.branch_means()
            _, ll = forward_backward(
                gaussian_log_emissions(xs, a, b, trial.sigma), p_flip)
            if ll > best_ll:
                best_ll, best_phi = float(ll), float(phi)
        step = (hi - lo) / max(n_grid - 1, 1)
        lo, hi = best_phi - step, best_phi + step
        n_grid = 9
    return float(np.mod(best_phi, ramp_period)), best_ll


def _cleanest_window_model(iq, dt, period, n_windows, kwargs):
    """Emission model whose profile comes from the single best window.

    A whole-trace profile is blended across any offset-charge jump; a per-window
    profile is clean for every window the jump misses. Blending always *costs*
    splitting contrast -- it averages two misaligned curves -- so picking the
    highest-contrast window reliably picks a clean one.
    """
    n = iq.size
    width = n // n_windows
    best = None
    for k in range(n_windows):
        model = learn_emission_model(
            iq, dt, fold_period=period,
            profile_slice=(k * width, (k + 1) * width), **kwargs)
        contrast = float(model.diagnostics.get("max_contrast", 0.0))
        if best is None or contrast > best[0]:
            best = (contrast, model)
    return best[1]


def reconstruct_parity_flips_ramped(
    iq: np.ndarray,
    sample_rate: float,
    t0: float = 0.0,
    p_flip_init: float = 1e-3,
    n_rate_iterations: int = 4,
    min_confidence: float = 0.0,
    emission: EmissionModel | None = None,
    segment_charge_jumps: bool = True,
    model_ramp_resets: bool = True,
    ramp_period: float | None = None,
    n_profile_windows: int = 6,
    n_segment_iterations: int = 3,
    min_detectability: float = 70.0,
    min_contrast: float = 1.0,
    decoder: str = "viterbi",
    **emission_kwargs,
) -> ReconstructionResult:
    """Recover parity-flip times from a complex readout trace, blind.

    No device or resonator parameters are used: the discrimination axis, the
    noise level, the two-branch splitting profile, the offset-charge ramp
    period and its reset schedule are all learned from ``iq`` itself. The
    branch sequence is then decoded with a two-state HMM and flips are read off
    the transitions.

    The pipeline is

    1. learn the projected two-branch emission model (:mod:`.emission`);
    2. locate offset-charge jumps and re-learn the model with the ramp phase
       re-aligned on each side (:mod:`.segment`);
    3. a short first-pass decode, whose event list is dominated by ramp-reset
       artefacts;
    4. find the reset comb in that list and install it in the model's sign
       schedule (:mod:`.ramp`), keeping it only if it improves the likelihood;
    5. final decode, and read off the transitions (:mod:`.events`).

    Parameters
    ----------
    iq : complex array
        Readout trace, ``I + 1j*Q``, on a uniform time grid.
    sample_rate : float
        Sampling rate [Hz].
    t0 : float
        Time of the first sample [s]; flip times are reported on this axis.
    p_flip_init : float
        Starting per-sample flip probability. The rate is re-estimated from the
        decoded path, so this only has to be within an order of magnitude.
    n_rate_iterations : int
        Hard-EM refinements of the flip rate (decode, count, re-decode).
    min_confidence : float
        Drop flips whose posterior swing is below this. The default keeps
        everything; raising it trades recall for precision, which mostly
        matters for flips landing in a parity-blind window.
    emission : EmissionModel, optional
        Pre-learned emission model; learned from ``iq`` when omitted. Supplying
        one skips steps 1-2.
    segment_charge_jumps : bool
        Re-align the ramp phase across offset-charge jumps. Jumps are treated
        purely as nuisances -- their amplitudes are not reported.
    ramp_period : float, optional
        The sawtooth period [s], when it is known from the hardware (i.e.
        ``1 / ramp_frequency``). The *phase* need not be known -- it is found
        by a likelihood scan. Supplying this replaces the blind comb search
        with something much more robust, and additionally re-locks the fold
        period to ``ramp_period / m`` for integer ``m``, which is the signal
        generator's precision rather than the fit's.

        This only matters at low contrast. Where the blind path already works
        (contrast of about 1.3 and above) it changes nothing; its value is that
        it moves the noise breakdown from a contrast near 0.9 down to about
        0.5, because comb detection -- not decoding -- is what fails first.
        If the supplied period is inconsistent with the trace the value is
        ignored and the blind search runs, with ``ramp_period_rejected`` set in
        the diagnostics.
    model_ramp_resets : bool
        Fold the sawtooth ramp-reset schedule into the emission model's sign
        schedule. A reset re-labels the branches without changing anything
        observable, so a decoder that does not know about it reports a flip at
        every reset -- 500 per second against a background tunnelling rate of
        order 10 Hz on the reference scenario.
    n_profile_windows : int
        Number of windows tried when picking a jump-free splitting profile to
        drive the change-point search (see :func:`_cleanest_window_model`).
    decoder : {"viterbi", "posterior"}
        Which decoding rule turns the HMM output into events; see
        :func:`~qpd.reconstruction.hmm.decoded_path`. Applies to the final
        extraction only -- the reset-comb search stays on Viterbi so the
        nuisance model a trace gets is independent of this choice.
    min_detectability, min_contrast : float
        Thresholds behind :attr:`ReconstructionResult.degenerate`. The trace is
        called unusable only when the median contrast falls below
        ``min_contrast`` *and* the contrast integrated over a dwell falls below
        ``min_detectability`` -- see the note there.
    n_segment_iterations : int
        Alternations between locating the jumps and polishing the fold period.
        One pass suffices for a fast ramp; a slow ramp has too few fold cycles
        for the windowed period estimate to shrug a jump off on its own and
        needs the second pass.

    Returns
    -------
    ReconstructionResult
    """
    iq = validate_trace(iq, sample_rate)
    dt = 1.0 / float(sample_rate)
    n = iq.size
    t = np.arange(n) * dt
    diag: dict = {}

    # A caller-supplied fold period must not also travel inside
    # emission_kwargs, or it collides with the explicit fold_period= that the
    # window-probe and re-learn steps pass.
    known_fold = emission_kwargs.pop("fold_period", None)

    jump_times = np.empty(0, dtype=float)
    if emission is None:
        emission = learn_emission_model(iq, dt, fold_period=known_fold,
                                        **emission_kwargs)

        if segment_charge_jumps and emission.fold_period is not None:
            # Period and jump offsets are coupled: a phase step biases the fold
            # objective, and a biased period drifts the phase in a way the
            # change-point search reads as a step. Alternate between them.
            x_proj = emission.project(iq)
            period = emission.fold_period
            phase_offset = np.zeros(n)
            jump_idx = np.empty(0, dtype=int)
            for _ in range(max(1, int(n_segment_iterations))):
                probe = _cleanest_window_model(
                    iq, dt, period, n_profile_windows, emission_kwargs)
                offs, idx = segment_and_realign(iq, dt, probe)
                if idx.size == 0:
                    break
                phase_offset, jump_idx = offs, idx
                # With the phase steps removed the whole trace folds coherently
                # again, so the period can be polished on the full duration.
                span = max(3e-5, 3.0 * period / (n * dt))
                new_period = refine_fold_period(
                    x_proj, dt, period, phase_offset, rel_span=span)
                converged = abs(new_period - period) < 1e-9 * period
                period = new_period
                if converged:
                    break
            if jump_idx.size:
                # Re-learn with the phase re-aligned: without the offsets,
                # samples from either side of a jump land in the same phase bin
                # and the learned splitting profile is the average of two
                # misaligned curves.
                emission = learn_emission_model(
                    iq, dt, phase_offset=phase_offset, fold_period=period,
                    **emission_kwargs)
                diag["segment_offsets"] = np.unique(phase_offset)
                diag["fold_period_polished"] = period
            jump_times = t0 + jump_idx * dt

    x = emission.project(iq)
    mu_a, mu_b = emission.branch_means()

    comb: ResetComb | None = None
    final: tuple | None = None  # (result, p, history) reused from the trial

    if (ramp_period is not None and model_ramp_resets
            and emission.fold_period is not None):
        ratio = float(ramp_period) / emission.fold_period
        m_fold = int(round(ratio))
        residual = abs(ratio - m_fold)
        diag["ramp_period_supplied"] = float(ramp_period)
        diag["ramp_period_fold_ratio"] = ratio
        if m_fold >= 2 and residual <= 0.15:
            # Lock the fold period to the hardware value: T/m is exact to the
            # generator's accuracy, better than the spectral fit, and it
            # removes the harmonic ambiguity the blind search has to resolve.
            locked = float(ramp_period) / m_fold
            # Keep the pre-lock model: if the supplied period fails
            # verification below we must hand the blind search the emission it
            # would have had, not one re-learnt at a wrong fold period.
            emission_before, x_before = emission, x
            if abs(locked - emission.fold_period) > 1e-12 * locked:
                emission = learn_emission_model(
                    iq, dt, phase_offset=emission.phase_offset,
                    fold_period=locked, **emission_kwargs)
                x = emission.project(iq)
            phi, _ = _scan_reset_phase(x, emission, t, float(ramp_period),
                                       p_flip_init)

            # Commensurability is NOT a test of the supplied period: almost any
            # value sits near some integer multiple of the fold period, so a
            # wrong one passes it and then silently produces garbage (measured:
            # T*1.37 accepted as n_fold=15, F1 0.002). The period has to earn
            # its place on evidence, by the same two tests the blind comb
            # faces -- it must improve the likelihood over using no comb at
            # all, and it must leave no periodic residue in the output.
            trial = emission.with_reset(t, float(ramp_period), phi)
            ta, tb = trial.branch_means()
            res_with, p_with, hist_with = decode_with_rate(
                x, ta, tb, trial.sigma, p_flip_init, n_rate_iterations)
            base, _, _ = decode_with_rate(
                x, mu_a, mu_b, emission.sigma, p_flip_init, 2)
            times_with, _ = extract_flips(res_with.path, res_with.posterior, dt)
            residue = (find_reset_comb(times_with, emission.fold_period,
                                       n * dt, dt)
                       if times_with.size >= 16 else None)
            better = res_with.log_likelihood > base.log_likelihood
            if better and residue is None:
                comb = ResetComb(period=float(ramp_period), phase=phi,
                                 n_fold=m_fold, excess=float("nan"),
                                 occupancy=float("nan"))
                emission = trial
                mu_a, mu_b = ta, tb
                final = (res_with, p_with, hist_with)
                diag.update({"fold_period_locked": locked,
                             "reset_phase_scanned": phi})
            else:
                diag["ramp_period_rejected"] = (
                    "supplied ramp_period did not survive verification "
                    f"(likelihood improved: {better}; periodic residue in the "
                    f"output: {residue is not None}); using the blind search")
                emission, x = emission_before, x_before
                mu_a, mu_b = emission.branch_means()
                ramp_period = None
        else:
            diag["ramp_period_rejected"] = (
                f"ramp_period/fold_period = {ratio:.3f} is not close to an "
                f"integer (nearest {m_fold}, residual {residual:.3f}); "
                "falling back to the blind comb search")
            ramp_period = None

    if comb is None and model_ramp_resets and emission.fold_period is not None:
        first, p_first, _ = decode_with_rate(
            x, mu_a, mu_b, emission.sigma, p_flip_init, 2)
        first_times, _ = extract_flips(first.path, first.posterior, dt)
        diag["first_pass_events"] = int(first_times.size)
        diag["first_pass_log_likelihood"] = first.log_likelihood
        candidate = find_reset_comb(first_times, emission.fold_period,
                                    n * dt, dt)
        if candidate is not None:
            trial = emission.with_reset(t, candidate.period, candidate.phase)
            ta, tb = trial.branch_means()
            res_t, p_t, hist_t = decode_with_rate(
                x, ta, tb, trial.sigma, p_flip_init, n_rate_iterations)
            # A reset schedule that is real explains away a large, strictly
            # periodic block of transitions, so it buys a lot of likelihood.
            # Keeping the gate means a trace with no such ramp is untouched.
            diag["reset_log_likelihood"] = res_t.log_likelihood
            if res_t.log_likelihood > first.log_likelihood:
                comb = candidate
                emission = trial
                mu_a, mu_b = ta, tb
                # Nothing downstream changes the emissions, so this decode *is*
                # the final one.
                final = (res_t, p_t, hist_t)
            else:
                diag["reset_rejected"] = True

    if final is None:
        final = decode_with_rate(
            x, mu_a, mu_b, emission.sigma, p_flip_init, n_rate_iterations)
    res, p, history = final

    # Only the FINAL event extraction honours `decoder`. The first-pass decode
    # and the comb trials above are nuisance detection, not output, and are
    # deliberately left on Viterbi so the reset schedule a trace gets does not
    # depend on which decoding rule the caller asked for.
    path = decoded_path(res, decoder)
    times, conf = extract_flips(path, res.posterior, dt, t0=t0)
    if min_confidence > 0:
        keep = conf >= min_confidence
        times, conf = times[keep], conf[keep]

    # Self-check: a periodic comb surviving in the *output* means the reset
    # schedule is wrong (typically locked to a multiple of the true period).
    residual = None
    if emission.fold_period is not None and times.size >= 16:
        residual = find_reset_comb(times - t0, emission.fold_period,
                                   n * dt, dt)
    starved = float(emission.diagnostics.get("phase_bin_counts_min", 1.0)) <= 0
    median_contrast = (float(np.median(emission.contrast))
                       if emission.separation.size else float("nan"))
    detectability = median_contrast * np.sqrt(1.0 / max(p, 1e-12))
    # Same two-condition test as the fixed-bias path: the branches must be
    # separated either per sample (contrast) or after integrating over a dwell
    # (detectability). Requiring both to fail is what stops a genuinely fast
    # telegraph being mislabelled, while still catching a trace that is simply
    # too noisy -- which the residual-comb test alone does not see, since a
    # decoder segmenting pure noise produces no periodic structure at all.
    too_noisy = bool(detectability < float(min_detectability)
                     and median_contrast < float(min_contrast))
    degenerate = bool(residual is not None or starved or too_noisy)

    diag.update({
        "residual_comb": residual,
        "decoder": decoder,
        "viterbi_path": res.path,
        "degenerate": degenerate,
        "detectability": float(detectability),
        "median_contrast": median_contrast,
        "starved_phase_grid": starved,
        "dt": dt,
        "rate_hz": p / dt,
        "p_flip_history": history,
        "n_flips": int(times.size),
        "log_likelihood": res.log_likelihood,
        "emission": emission.diagnostics,
        "reset_comb": comb,
    })
    return ReconstructionResult(
        flip_times=times,
        confidence=conf,
        posterior=res.posterior,
        branch=path,
        emission=emission,
        p_flip=p,
        reset_comb=comb,
        charge_jump_times=jump_times,
        diagnostics=diag,
    )
