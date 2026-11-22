"""Timing and preservation diagnostics for the LED campaign."""
import numpy as np

from .led_generalization_score import match_events, timing_settings


def strict_event_count(truth, predicted, tolerance):
    """Maximum count of one-to-one matches within an inclusive time tolerance.

    Both event lists must be sorted. If the earliest events are compatible,
    matching them preserves an optimal assignment. Otherwise the earlier
    event cannot match anything in the later list.
    """
    truth, predicted = (np.asarray(a, dtype=float) for a in (truth, predicted))
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("Expected a finite nonnegative tolerance")
    for times in (truth, predicted):
        if times.ndim != 1 or not np.all(np.isfinite(times)) or np.any(np.diff(times) < 0):
            raise ValueError("Expected finite, sorted, one-dimensional event times")
    i = j = count = 0
    while i < len(truth) and j < len(predicted):
        delta = predicted[j] - truth[i]
        if delta < -tolerance:
            j += 1
        elif delta > tolerance:
            i += 1
        else:
            count += 1
            i += 1
            j += 1
    return count


def c1_truth_states(truth, predicted, fs):
    """Per-truth status under frozen whole-record C1: 0 unmatched, 1 late, 2 precise."""
    ti, _, delta = match_events(truth, predicted, fs)
    states = np.zeros(len(truth), dtype=np.int8)
    states[ti] = 1 + (np.abs(delta) <= timing_settings(fs)[0]).astype(np.int8)
    return states


def preservation_counts(before, after):
    """3-by-3 paired truth transition counts, rows baseline and columns candidate."""
    if len(before) != len(after):
        raise ValueError("Both status arrays must describe the same truth events")
    return np.bincount(3*np.asarray(before) + np.asarray(after), minlength=9).reshape(3, 3)
