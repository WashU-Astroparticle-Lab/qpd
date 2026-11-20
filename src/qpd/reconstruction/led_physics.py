"""Calibrated charge/parity smoothing and timestamp selection for LED records."""
from dataclasses import dataclass
from operator import index
from typing import Optional

import numpy as np
from scipy.stats import exponnorm

from .benchmark import as_complex_trace


@dataclass
class LEDPhysicsResult:
    flip_times: np.ndarray
    # [gap, subinterval], in expected physical flips (not toggle probabilities).
    expected_counts: np.ndarray
    # [sample, calibration column], when requested.
    charge_posterior: Optional[np.ndarray]


def _expected_counts(times, rises, on, quiet, mean, mu, sigma, tau, subdivisions):
    step = (times[1] - times[0]) / subdivisions
    boundaries = np.arange((len(times) - 1) * subdivisions + 1) * step
    counts = np.full(len(boundaries) - 1, quiet * step)
    if mean:
        law = exponnorm(tau / sigma, loc=mu, scale=sigma)
        zero = law.cdf(0.)
        survival = law.sf(0.)
        # The omitted conditional tail is at most 1e-14 per LED event.
        horizon = on + law.isf(1e-14 * survival)
        nodes, weights = np.polynomial.legendre.leggauss(8)
        for rise in rises:
            left = max(0, boundaries.searchsorted(rise) - 1)
            right = boundaries.searchsorted(rise + horizon) + 1
            lag = np.maximum(0., boundaries[left:right] - rise)
            width = np.minimum(lag, on)
            cdf = np.zeros_like(lag)
            for node, weight in zip(nodes, weights):
                cdf += weight * (law.cdf(lag - width * (1 + node) / 2) - zero)
            cdf *= width / (2 * on * survival)
            counts[left:left + len(lag) - 1] += mean * np.maximum(0., np.diff(cdf))
    return counts.reshape(-1, subdivisions)


def _align(iq, branches, sigma):
    """Fit one rigid rotation and translation from the measured point cloud."""
    cloud = branches.ravel()
    y = iq[::8]
    y_centered = y - y.mean()
    axis = .5 * np.angle(np.mean(y_centered**2) * np.conj(np.mean((cloud-cloud.mean())**2)))
    best_score = -np.inf
    for initial in (axis, axis + np.pi):
        rotation = np.exp(1j * initial)
        shift = y.mean() - rotation * cloud.mean()
        for iteration in range(7):
            logits = -.5 * (abs(y[:, None] - shift - rotation * cloud) / sigma)**2
            peak = logits.max(axis=1)
            weights = np.exp(logits - peak[:, None])
            sums = weights.sum(axis=1)
            if iteration == 6:
                score = (peak + np.log(sums)).sum()
                if score > best_score:
                    best_score = score
                    aligned = (iq - shift) / rotation
                break
            centers = (weights @ cloud) / sums
            rotation = np.exp(1j * np.angle(np.vdot(centers-centers.mean(), y_centered)))
            shift = y.mean() - rotation * centers.mean()
    return aligned


def _smooth(iq, branches, sigma, lambdas, redraw, return_charge):
    # Log space preserves unlikely paths across very long quiet runs or a
    # zero-rate gap. Both transition kernels are symmetric and commute.
    switch = -.5 * np.expm1(-2 * lambdas)
    with np.errstate(divide='ignore'):
        log_switch = np.log(switch)
        log_stay = np.log1p(-switch)
        log_reset = np.log(redraw) - np.log(branches.shape[1])
        log_hold = np.log1p(-redraw)

    def mix_charge(state, k):
        total = np.logaddexp.reduce(state, axis=1, keepdims=True)
        return np.logaddexp(state + log_hold[k], total + log_reset[k])

    def emission(k):
        return -.5 * (abs(iq[k] - branches) / sigma)**2

    forward = np.empty((len(iq), *branches.shape))
    initial = emission(0)  # uniform prior over parity and charge
    forward[0] = initial - initial.max()
    for k in range(1, len(iq)):
        prev = forward[k-1]
        state = np.logaddexp(prev + log_stay[k-1], prev[::-1] + log_switch[k-1])
        state = mix_charge(state, k) + emission(k)
        forward[k] = state - state.max()

    charge = np.empty((len(iq), branches.shape[1])) if return_charge else None
    backward = np.zeros_like(branches.real)
    odd = np.empty_like(lambdas)
    for k in range(len(iq)-1, -1, -1):
        if return_charge:
            log_joint = forward[k] + backward
            marginal = np.exp(log_joint - log_joint.max()).sum(axis=0)
            charge[k] = marginal / marginal.sum()
        if k == 0:
            break
        after = mix_charge(backward + emission(k), k)
        prev = forward[k-1]
        cross = log_switch[k-1] + np.logaddexp.reduce((prev + after[::-1]).ravel())
        stay = log_stay[k-1] + np.logaddexp.reduce((prev + after).ravel())
        odd[k-1] = np.exp(cross - np.logaddexp(cross, stay))
        backward = np.logaddexp(after + log_stay[k-1], after[::-1] + log_switch[k-1])
        backward -= backward.max()
    return odd, charge


def _timestamps(counts, step, tolerance, budget, threshold):
    mass = counts.ravel()
    width = min(len(mass), max(1, round(2 * tolerance / step)))
    cumulative = np.r_[0., mass.cumsum()]
    gains = cumulative[width:] - cumulative[:-width]
    blocked = np.zeros(len(gains), dtype=bool)
    calls = []
    for start in np.argsort(-gains, kind='stable'):
        if len(calls) == budget or gains[start] < threshold:
            break
        if not blocked[start]:
            calls.append((start + width / 2) * step)
            blocked[max(0, start-width+1):start+width] = True
    return np.sort(calls)


def reconstruct_parity_flips_led_physics(
    iq, sample_rate, first_pulse_time_s, *, on_time_s, off_time_s,
    branch_iq, noise_sigma, call_budget=None, min_expected_count=None,
    tolerance_s=None, quiescent_rate_hz=20., mean_led_flips=1.,
    burst_mu_s=.0002, burst_sigma_s=.0004, burst_tau_s=.0008,
    charge_jump_rate_hz=1/.6, align_iq=False, return_charge=False,
):
    """Infer continuous-time physical flips from sampled complex I/Q.

    Times are seconds relative to sample zero. The first physical LED rise is
    first_pulse_time_s (no earlier pulses); ON plus OFF is the period.
    branch_iq has shape (2, n_charge): even and odd complex means at the same
    equally weighted charge bins. noise_sigma is the single-sample standard
    deviation per quadrature. Set align_iq to estimate one rigid I/Q transform
    when calibration and waveform use different orientations or offsets.

    LED counts are Poisson. Excitation is uniform during ON with nonnegative-
    conditioned EMG delay. Each rise redraws charge; other redraws are Poisson.
    Select exactly one operating policy: call_budget (maximum calls per record)
    or min_expected_count (positive expected mass per selected window).
    Greedy disjoint windows can exhaust available space before spending a budget.
    Tolerance defaults to max(50 us, 2/sample_rate). Window widths are rounded
    to the nearest eighth of a sample, with at least one subinterval.
    charge_posterior columns correspond to the supplied calibration columns.
    """
    iq = as_complex_trace(iq)
    branches = np.asarray(branch_iq, dtype=complex)
    if len(iq) < 2 or not np.isfinite(iq).all():
        raise ValueError("expected at least two finite I/Q samples")
    if branches.ndim != 2 or branches.shape[0] != 2 or not branches.shape[1] or not np.isfinite(branches).all():
        raise ValueError("branch_iq must contain finite paired means of shape (2, n_charge)")
    positive = (sample_rate, on_time_s, off_time_s, noise_sigma, burst_sigma_s, burst_tau_s)
    nonnegative = (first_pulse_time_s, mean_led_flips, quiescent_rate_hz, burst_mu_s, charge_jump_rate_hz)
    if not all(np.isfinite(x) and x > 0 for x in positive):
        raise ValueError("sample rate, ON/OFF durations, noise and delay scales must be positive and finite")
    if not all(np.isfinite(x) and x >= 0 for x in nonnegative):
        raise ValueError("onset, rates, burst location and LED mean must be nonnegative and finite")
    tolerance = max(50e-6, 2/sample_rate) if tolerance_s is None else tolerance_s
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance_s must be positive and finite")
    if (call_budget is None) == (min_expected_count is None):
        raise ValueError("specify exactly one of call_budget or min_expected_count")
    if call_budget is not None and index(call_budget) < 0:
        raise ValueError("call_budget must be a nonnegative integer")
    if min_expected_count is not None and (not np.isfinite(min_expected_count) or min_expected_count <= 0):
        raise ValueError("min_expected_count must be positive and finite")

    times = np.arange(len(iq)) * (1/sample_rate)
    period = on_time_s + off_time_s
    n_rises = max(0, int(np.floor((times[-1] - first_pulse_time_s)/period)) + 1)
    rises = first_pulse_time_s + np.arange(n_rises) * period
    rises = rises[rises <= times[-1]]
    subdivisions = 8
    counts = _expected_counts(times, rises, on_time_s, quiescent_rate_hz,
                              mean_led_flips, burst_mu_s, burst_sigma_s, burst_tau_s, subdivisions)
    redraw = np.full(len(iq), -np.expm1(-charge_jump_rate_hz/sample_rate))
    redraw[times.searchsorted(rises)] = 1.
    if align_iq:
        iq = _align(iq, branches, noise_sigma)
    lambdas = counts.sum(axis=1)
    odd, charge = _smooth(iq, branches, noise_sigma, lambdas, redraw, return_charge)
    # Conditional Poisson means for even/odd counts in a gap. Inside a gap
    # the conditional event-time density stays proportional to the prior rate.
    tanh = np.tanh(lambdas)
    scale = (1-odd)*tanh + np.divide(odd, tanh, out=np.zeros_like(odd), where=tanh > 0)
    counts *= scale[:, None]
    calls = _timestamps(counts, 1/(subdivisions*sample_rate), tolerance,
                        call_budget, min_expected_count or 0.)
    return LEDPhysicsResult(calls, counts, charge)
