"""End-to-end blind reconstruction of parity-flip timing from an I/Q trace."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .emission import EmissionModel, learn_emission_model, refine_fold_period
from .events import extract_flips
from .hmm import HMMResult, forward_backward, gaussian_log_emissions, viterbi
from .ramp import ResetComb, find_reset_comb
from .segment import segment_and_realign

__all__ = ["ReconstructionResult", "reconstruct_parity_flips"]


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
    def contrast(self) -> np.ndarray:
        """Per-sample branch separation in units of the noise std."""
        return (self.emission.contrast if self.emission is not None
                else np.empty(0))


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


def _decode_with_rate(x, mu_a, mu_b, sigma, p0, n_iter):
    """Decode, re-estimating the flip rate from the decoded path (hard EM).

    The rate only sets the decoder's prior, so a handful of hard-EM passes is
    plenty; a full Baum-Welch update would cost a second sweep per iteration for
    no measurable gain here. The emissions are built once, and the iterations
    run forward-backward only -- the Viterbi path is needed just once, at the
    end, to enumerate the events.
    """
    n = x.size
    log_emit = gaussian_log_emissions(x, mu_a, mu_b, sigma)
    p = float(p0)
    history = [p]
    for _ in range(max(1, int(n_iter))):
        post, _ = forward_backward(log_emit, p)
        # Thresholded posterior stands in for the Viterbi path here: it is a
        # transition count, and the two agree except inside parity-blind
        # windows where the posterior sits near 1/2 either way.
        crossings = np.count_nonzero(np.diff(post > 0.5))
        p_new = max(crossings / max(n - 1, 1), 1e-9)
        history.append(p_new)
        converged = abs(p_new - p) < 1e-3 * p
        p = p_new
        if converged:
            break
    post, loglik = forward_backward(log_emit, p)
    path = viterbi(log_emit, p)
    res = HMMResult(posterior=post, path=path, log_likelihood=loglik, p_flip=p,
                    diagnostics={"p_flip_history": history})
    return res, p, history


def reconstruct_parity_flips(
    iq: np.ndarray,
    sample_rate: float,
    t0: float = 0.0,
    p_flip_init: float = 1e-3,
    n_rate_iterations: int = 4,
    min_confidence: float = 0.0,
    emission: EmissionModel | None = None,
    segment_charge_jumps: bool = True,
    model_ramp_resets: bool = True,
    n_profile_windows: int = 6,
    n_segment_iterations: int = 3,
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
    model_ramp_resets : bool
        Fold the sawtooth ramp-reset schedule into the emission model's sign
        schedule. A reset re-labels the branches without changing anything
        observable, so a decoder that does not know about it reports a flip at
        every reset -- 500 per second against a background tunnelling rate of
        order 10 Hz on the reference scenario.
    n_profile_windows : int
        Number of windows tried when picking a jump-free splitting profile to
        drive the change-point search (see :func:`_cleanest_window_model`).
    n_segment_iterations : int
        Alternations between locating the jumps and polishing the fold period.
        One pass suffices for a fast ramp; a slow ramp has too few fold cycles
        for the windowed period estimate to shrug a jump off on its own and
        needs the second pass.

    Returns
    -------
    ReconstructionResult
    """
    iq = np.asarray(iq)
    if iq.ndim != 1 or iq.size < 3:
        raise ValueError("iq must be a 1-D trace with at least 3 samples")
    dt = 1.0 / float(sample_rate)
    n = iq.size
    t = np.arange(n) * dt
    diag: dict = {}

    jump_times = np.empty(0, dtype=float)
    if emission is None:
        emission = learn_emission_model(iq, dt, **emission_kwargs)

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
    if model_ramp_resets and emission.fold_period is not None:
        first, p_first, _ = _decode_with_rate(
            x, mu_a, mu_b, emission.sigma, p_flip_init, 2)
        first_times, _ = extract_flips(first.path, first.posterior, dt)
        diag["first_pass_events"] = int(first_times.size)
        diag["first_pass_log_likelihood"] = first.log_likelihood
        candidate = find_reset_comb(first_times, emission.fold_period,
                                    n * dt, dt)
        if candidate is not None:
            trial = emission.with_reset(t, candidate.period, candidate.phase)
            ta, tb = trial.branch_means()
            res_t, p_t, hist_t = _decode_with_rate(
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
        final = _decode_with_rate(
            x, mu_a, mu_b, emission.sigma, p_flip_init, n_rate_iterations)
    res, p, history = final

    times, conf = extract_flips(res.path, res.posterior, dt, t0=t0)
    if min_confidence > 0:
        keep = conf >= min_confidence
        times, conf = times[keep], conf[keep]

    diag.update({
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
        branch=res.path,
        emission=emission,
        p_flip=p,
        reset_comb=comb,
        charge_jump_times=jump_times,
        diagnostics=diag,
    )
