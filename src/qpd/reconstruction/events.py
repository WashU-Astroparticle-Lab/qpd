"""Turn a decoded branch sequence into timed parity-flip events."""

from __future__ import annotations

import numpy as np

__all__ = ["extract_flips", "flip_confidence"]


def extract_flips(
    path: np.ndarray,
    posterior: np.ndarray,
    dt: float,
    t0: float = 0.0,
    refine: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Flip times and per-flip confidences from a decoded branch sequence.

    A flip sits between the last sample of one branch and the first of the
    next, so the default estimate is the midpoint of that interval. When
    ``refine`` is set the time is instead read off where the posterior odds
    cross even, linearly interpolated between the bracketing samples: with a
    10 us sample period against a sub-millisecond timing tolerance this is a
    small correction, but it is free and unbiased.

    Confidence is the posterior swing across the transition, in [0, 1]. It
    collapses toward zero for flips inside a parity-blind window, where the
    data genuinely carries no information about the branch.
    """
    path = np.asarray(path)
    posterior = np.asarray(posterior, dtype=float)
    idx = np.flatnonzero(np.diff(path) != 0) + 1
    if idx.size == 0:
        return np.empty(0, dtype=float), np.empty(0, dtype=float)

    times = t0 + (idx - 0.5) * dt
    if refine:
        p_before = posterior[idx - 1]
        p_after = posterior[idx]
        # Linear interpolation of the 0.5 level between the two samples; only
        # meaningful when the pair actually brackets it.
        span = p_after - p_before
        brackets = np.abs(span) > 1e-12
        frac = np.where(brackets,
                        (0.5 - p_before) / np.where(brackets, span, 1.0), 0.5)
        frac = np.clip(frac, 0.0, 1.0)
        times = t0 + (idx - 1 + frac) * dt

    return times, flip_confidence(posterior, idx)


def flip_confidence(posterior: np.ndarray, idx: np.ndarray) -> np.ndarray:
    """Posterior swing across each transition, clipped to [0, 1]."""
    return np.clip(np.abs(posterior[idx] - posterior[idx - 1]), 0.0, 1.0)
