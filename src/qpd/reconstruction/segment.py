"""Ramp-phase re-alignment across offset-charge jumps.

A parity flip and an offset-charge jump look different in this readout. A flip
moves the signal *between* two branch curves that themselves stay put. A jump
of ``delta`` in n_g slides both curves along the ramp phase by ``2*delta``, so
the learned curves are wrong from the jump onwards -- transient versus
persistent.

That asymmetry is what lets jumps be handled as nuisances: find where the
ramp phase steps, re-fit one phase offset per segment, and emit no flip at the
boundary. Jump amplitudes are never reported; they are not a reconstruction
target here.

Ignoring jumps is not an option even when they are rare. The phase offset
enters the *sign* schedule as well as the splitting magnitude, so a jump that
lands near ``delta = 0.25`` (a half-cycle shift) inverts the model's idea of
which branch is which and the decoder degenerates into flipping constantly.
On the reference scenario a single such jump in a 5 s trace drops the flip F1
from 0.98 to 0.21.

One case is irreducible: ``delta = +/-0.5`` (mod 1) maps each branch exactly
onto the other, which is indistinguishable from a parity flip in any single
trace. The same applies to the sawtooth reset, which is handled separately and
statistically in :mod:`.ramp`.
"""

from __future__ import annotations

import numpy as np

from .emission import EmissionModel

__all__ = ["segment_and_realign", "misfit_statistic"]


def misfit_statistic(x: np.ndarray, mu_a: np.ndarray, mu_b: np.ndarray,
                     sigma: float) -> np.ndarray:
    """Squared distance to the nearer branch, in units of the noise variance.

    Close to 1 (in fact somewhat below, being the smaller of two draws) when
    the model is right, and inflated wherever the ramp phase is wrong. Note it
    depends only on the *pair* of branch means, not on which is which, so it is
    blind to the sign schedule -- and therefore usable before the ramp resets
    have been found.
    """
    da = (x - mu_a) ** 2
    db = (x - mu_b) ** 2
    return np.minimum(da, db) / (sigma * sigma)


def _block_offsets(x: np.ndarray, t: np.ndarray, model: EmissionModel,
                   block: int, n_offsets: int) -> np.ndarray:
    """Best ramp-phase offset per contiguous block of samples."""
    n = x.size
    starts = np.arange(0, n, block)
    grid = np.arange(n_offsets) / n_offsets
    best = np.empty(starts.size)
    for b, lo in enumerate(starts):
        hi = min(lo + block, n)
        xs, ts = x[lo:hi], t[lo:hi]
        cost = np.empty(n_offsets)
        for j, off in enumerate(grid):
            c, h = model.evaluate(ts, float(off))
            cost[j] = np.mean(
                misfit_statistic(xs, c + 0.5 * h, c - 0.5 * h, model.sigma)
            )
        best[b] = grid[int(np.argmin(cost))]
    return best


def _refine_boundary(x: np.ndarray, t: np.ndarray, model: EmissionModel,
                     lo: int, hi: int, off_before: float,
                     off_after: float) -> int:
    """Locate the change-point inside ``[lo, hi)`` by cumulative misfit."""
    xs, ts = x[lo:hi], t[lo:hi]
    ca, ha = model.evaluate(ts, off_before)
    cb, hb = model.evaluate(ts, off_after)
    ma = misfit_statistic(xs, ca + 0.5 * ha, ca - 0.5 * ha, model.sigma)
    mb = misfit_statistic(xs, cb + 0.5 * hb, cb - 0.5 * hb, model.sigma)
    # Total cost of switching models at each candidate index.
    total = np.concatenate(([0.0], np.cumsum(ma))) + (
        np.sum(mb) - np.concatenate(([0.0], np.cumsum(mb)))
    )
    return lo + int(np.argmin(total))


def _circular_median(v: np.ndarray) -> float:
    """Median of phases on the unit circle (period 1)."""
    if v.size == 1:
        return float(v[0])
    ang = 2 * np.pi * v
    ref = np.angle(np.mean(np.exp(1j * ang)))
    centred = (ang - ref + np.pi) % (2 * np.pi) - np.pi
    return float(np.mod((ref + np.median(centred)) / (2 * np.pi), 1.0))


def _circular_median_filter(v: np.ndarray, width: int) -> np.ndarray:
    """Running circular median, for suppressing single-block phase noise."""
    if width <= 1 or v.size <= width:
        return v
    half = width // 2
    out = np.empty_like(v)
    for i in range(v.size):
        out[i] = _circular_median(v[max(0, i - half):min(v.size, i + half + 1)])
    return out


def _prune_steps(steps: np.ndarray, n_blocks: int, min_gap: int) -> np.ndarray:
    """Keep only well-separated change-points, away from the trace ends.

    A real offset-charge jump is a rare, permanent event (~0.1 Hz), so genuine
    boundaries are far apart. Phase noise, by contrast, produces bursts of
    candidate steps that flip-flop between two offsets over a few blocks; that
    is what this drops. Without it, a trace with long parity dwells -- where each
    block sits on a single branch and the ramp phase is only weakly constrained
    -- can be split into a dozen bogus segments.
    """
    kept: list[int] = []
    for s in steps:
        if s < min_gap or s > n_blocks - 1 - min_gap:
            continue
        if kept and s - kept[-1] < min_gap:
            continue
        kept.append(int(s))
    return np.asarray(kept, dtype=int)


def segment_and_realign(
    iq: np.ndarray,
    dt: float,
    model: EmissionModel,
    min_block_samples: int = 1024,
    block_periods: float = 4.0,
    n_offsets: int = 48,
    min_offset_change: float = 0.08,
    median_blocks: int = 5,
    min_segment_blocks: int = 16,
) -> tuple[np.ndarray, np.ndarray]:
    """Find offset-charge jumps and the ramp-phase offset of each segment.

    Returns ``(phase_offset, boundaries)``: a per-sample offset array in
    fold-period units, and the sample indices where a new segment starts. When
    no jump is found the offset array is all zeros and ``boundaries`` is empty.

    The block length is floored at ``min_block_samples`` rather than being a
    fixed number of fold periods: at a 500 Hz ramp four fold periods is only
    ~70 samples, far too few to pin a phase at a per-sample contrast of order
    unity, and the resulting noise would manufacture spurious jumps.

    The search is deliberately conservative -- a running circular median over
    ``median_blocks`` blocks, then a minimum separation of ``min_segment_blocks``
    between accepted boundaries. A false jump is expensive (it installs a wrong
    ramp phase over a whole segment) while a missed one costs only the contrast
    lost to a blended profile, so the asymmetry is intentional. One consequence
    worth knowing: a jump within ``min_segment_blocks`` of either end of the
    trace is left uncorrected, since there is too little data on the short side
    to fit its offset reliably.
    """
    n = iq.size
    zeros = np.zeros(n)
    if model.fold_period is None or model.magnitude_profile is None:
        return zeros, np.empty(0, dtype=int)

    x = model.project(iq)
    t = np.arange(n) * dt
    block = max(int(round(block_periods * model.fold_period / dt)),
                int(min_block_samples))
    if block >= n:
        return zeros, np.empty(0, dtype=int)

    best = _block_offsets(x, t, model, block, n_offsets)
    best = _circular_median_filter(best, median_blocks)

    # A jump shows up as a step in the per-block offset. Compare on the
    # circle: the offset is a phase, so 0.98 and 0.02 are neighbours.
    d = np.abs((np.diff(best) + 0.5) % 1.0 - 0.5)
    steps = _prune_steps(np.flatnonzero(d > min_offset_change), best.size,
                         min_segment_blocks)
    if steps.size == 0:
        return zeros, np.empty(0, dtype=int)

    # Consolidate: circular median offset over each run of blocks between steps.
    edges = np.concatenate(([0], steps + 1, [best.size]))
    offsets = np.array([
        _circular_median(best[a:b]) for a, b in zip(edges[:-1], edges[1:])
    ])
    boundaries = []
    for s in steps:
        lo = int(s * block)
        hi = min(int((s + 2) * block), n)
        seg = int(np.searchsorted(steps, s, side="left"))
        boundaries.append(
            _refine_boundary(x, t, model, lo, hi,
                             float(offsets[seg]), float(offsets[seg + 1]))
        )
    boundaries = np.asarray(boundaries, dtype=int)

    phase_offset = np.empty(n)
    edges_s = np.concatenate(([0], boundaries, [n]))
    for s, (lo, hi) in enumerate(zip(edges_s[:-1], edges_s[1:])):
        phase_offset[lo:hi] = offsets[s]
    # Report offsets relative to the first segment, so a jump-free trace keeps
    # the zero it started with.
    return phase_offset - offsets[0], boundaries
