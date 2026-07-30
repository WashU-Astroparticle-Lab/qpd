"""Detection of the sawtooth ramp-reset schedule.

A sawtooth ``n_g`` ramp spanning an odd number of half Cooper pairs (5.5 for
the reference device) steps ``n_g`` by half a Cooper pair modulo one every
time it resets. That maps each parity branch exactly onto the other, so a
reset is *observationally indistinguishable* from a parity flip: no per-event
test can separate them.

What can be used is that resets are strictly periodic while tunnelling is
Poisson. At a 500 Hz ramp there are 500 resets per second against a
background tunnelling rate of order 10 Hz, so the resets dominate a first-pass
event list by more than an order of magnitude and their phase pile-up is
overwhelming.

Removing the detected events afterwards is the wrong move at this ramp rate:
one removal per 2 ms cycle, with a matching window comparable to the timing
tolerance, would delete a substantial fraction of the real flips too. Instead
the comb found here is installed in the emission model's sign schedule (see
:meth:`.EmissionModel.reset_sign`), so the decoder never reports those
transitions in the first place, and no real event is ever at risk.

Because the reset interval is an integer number of half Cooper pairs of ramp,
the search is restricted to integer multiples of the already-known fold
period, which makes it both cheap and sharp.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["ResetComb", "find_reset_comb"]


@dataclass
class ResetComb:
    """A strictly periodic ramp-reset schedule."""

    period: float  # [s] ramp period
    phase: float  # [s] reset time within one period
    n_fold: int  # ramp period in units of the fold period
    excess: float  # significance of the phase pile-up, in sigma
    occupancy: float  # detected events per ramp cycle


def _pileup(times: np.ndarray, period: float, tol: float):
    """Strongest phase pile-up at ``period``: (excess, phase, count)."""
    n_bins = max(int(round(period / tol)), 4)
    idx = np.minimum((np.mod(times, period) / period * n_bins).astype(np.intp),
                     n_bins - 1)
    counts = np.bincount(idx, minlength=n_bins)
    b = int(np.argmax(counts))
    expected = times.size / n_bins
    excess = (counts[b] - expected) / np.sqrt(max(expected, 1.0))
    return float(excess), (b + 0.5) / n_bins * period, int(counts[b])


def _refit(times: np.ndarray, period: float, phase: float,
           dt: float) -> tuple[float, float]:
    """Least-squares refit of period and phase to the comb members.

    The pile-up bin locates the period only to the bin width; fitting a line
    through the member times against their cycle index uses the full lever arm
    of the trace, which matters because the schedule has to stay aligned over
    thousands of cycles.
    """
    for window in (4.0 * dt, 2.0 * dt, 1.0 * dt):
        cycle = np.round((times - phase) / period)
        resid = times - phase - cycle * period
        sel = np.abs(resid) < window
        if int(sel.sum()) < 8:
            break
        design = np.vstack([np.ones(int(sel.sum())), cycle[sel]]).T
        sol, *_ = np.linalg.lstsq(design, times[sel], rcond=None)
        phase, period = float(sol[0]), float(sol[1])
    return period, phase


def find_reset_comb(
    times: np.ndarray,
    fold_period: float,
    duration: float,
    dt: float,
    n_fold_max: int = 64,
    min_excess: float = 8.0,
    min_occupancy: float = 0.35,
    max_occupancy: float = 1.4,
) -> ResetComb | None:
    """Find the ramp-reset comb hiding in a first-pass event list.

    Candidate periods are integer multiples ``k * fold_period``: the ramp spans
    a whole number of half Cooper pairs, so the reset interval is a whole
    number of fold periods. Each candidate is scored by the significance of
    its strongest phase pile-up.

    The occupancy window is what rejects multiples of the true period. Folding
    at ``2 T`` collects two resets per cycle and so scores an occupancy near 2,
    while the true period lands near 1 -- one reset per ramp.

    Returns None when nothing stands out, which is the correct answer for a
    static ``n_g``, for an even-span sawtooth (whose resets leave the branch
    labels alone), or for a triangle ramp with no reset at all.
    """
    times = np.sort(np.asarray(times, dtype=float))
    if times.size < 16 or fold_period <= 0:
        return None

    best: ResetComb | None = None
    for k in range(2, n_fold_max + 1):
        period = k * fold_period
        if period > duration / 8.0:
            break
        excess, phase, count = _pileup(times, period, 2.0 * dt)
        occupancy = count / max(duration / period, 1.0)
        if not (min_occupancy <= occupancy <= max_occupancy):
            continue
        if best is None or excess > best.excess:
            best = ResetComb(period=period, phase=phase, n_fold=k,
                             excess=excess, occupancy=occupancy)

    if best is None or best.excess < min_excess:
        return None

    period, phase = _refit(times, best.period, best.phase, dt)
    if period <= 0:
        return None
    return ResetComb(period=period, phase=float(np.mod(phase, period)),
                     n_fold=best.n_fold, excess=best.excess,
                     occupancy=best.occupancy)
