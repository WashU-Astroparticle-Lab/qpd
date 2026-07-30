"""Scoring a reconstruction against truth, resolved by time window.

The scalar figures of merit come from :mod:`qpd.mlebench.grade`, which is the
grader the datasets are scored with: predicted times are matched one-to-one to
truth times by minimum total timing error subject to a hard tolerance.

What this module adds is the *windowed* view needed to study a quasiparticle
burst. Efficiency and timing bias both depend on how many tunnelling events
land in the window: a burst crowds tens of flips into a few milliseconds, and
once two flips fall inside the same handful of samples the trace cannot
resolve them at all. :func:`window_report` measures that directly by counting
truth events per window rather than assuming a rate.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from qpd.mlebench.grade import match_times

__all__ = ["FlipScore", "score_flips", "window_report", "burst_report"]

DEFAULT_TOL = 0.5e-3  # [s] hard matching tolerance, as used by the grader
DEFAULT_SCALE = 0.25e-3  # [s] timing kernel scale for the soft score


@dataclass
class FlipScore:
    """Detection and timing performance on one set of flip times."""

    n_truth: int
    n_pred: int
    n_matched: int
    hard_f1: float
    soft_f1: float
    efficiency: float  # matched / truth  (recall)
    purity: float  # matched / predicted (precision)
    bias_s: float  # mean signed timing residual (pred - truth)
    rms_s: float  # rms timing residual
    dt: np.ndarray = field(default_factory=lambda: np.empty(0, float))

    def row(self) -> str:
        return (f"{self.n_truth:6d} {self.n_pred:6d} {self.n_matched:6d} "
                f"{self.efficiency:8.3f} {self.purity:7.3f} "
                f"{self.hard_f1:7.3f} {self.soft_f1:7.3f} "
                f"{self.bias_s * 1e6:+9.1f} {self.rms_s * 1e6:8.1f}")

    @staticmethod
    def header() -> str:
        return (f"{'truth':>6s} {'pred':>6s} {'match':>6s} {'effic':>8s} "
                f"{'purity':>7s} {'hardF1':>7s} {'softF1':>7s} "
                f"{'bias/us':>9s} {'rms/us':>8s}")


def score_flips(truth: np.ndarray, pred: np.ndarray, tol: float = DEFAULT_TOL,
                scale: float = DEFAULT_SCALE) -> FlipScore:
    """Match predicted to true flip times and summarise the result."""
    truth = np.asarray(truth, dtype=float)
    pred = np.asarray(pred, dtype=float)
    m = match_times(truth, pred, tol)
    n_m = int(m.dt.size)
    q = np.exp(-(m.dt / scale) ** 2) if n_m else np.zeros(0)
    soft_tp = float(q.sum())
    eff = n_m / m.n_truth if m.n_truth else np.nan
    pur = n_m / m.n_pred if m.n_pred else np.nan
    denom = m.n_truth + m.n_pred
    sp = soft_tp / m.n_pred if m.n_pred else 0.0
    sr = soft_tp / m.n_truth if m.n_truth else 0.0
    return FlipScore(
        n_truth=int(m.n_truth), n_pred=int(m.n_pred), n_matched=n_m,
        hard_f1=(2 * n_m / denom) if denom else np.nan,
        soft_f1=(2 * sp * sr / (sp + sr)) if (sp + sr) > 0 else 0.0,
        efficiency=eff, purity=pur,
        bias_s=float(np.mean(m.dt)) if n_m else np.nan,
        rms_s=float(np.sqrt(np.mean(m.dt ** 2))) if n_m else np.nan,
        dt=m.dt,
    )


def window_report(
    truth: np.ndarray,
    pred: np.ndarray,
    windows: np.ndarray,
    tol: float = DEFAULT_TOL,
    scale: float = DEFAULT_SCALE,
) -> list[tuple[tuple[float, float], FlipScore]]:
    """Score each ``(t_start, t_end)`` window separately.

    Matching is run inside each window after padding by ``tol``, so an event
    matched across the edge is not lost, and predictions are attributed to the
    window their truth partner belongs to.
    """
    out = []
    for lo, hi in np.asarray(windows, dtype=float).reshape(-1, 2):
        tsel = (truth >= lo) & (truth < hi)
        psel = (pred >= lo - tol) & (pred < hi + tol)
        out.append(((float(lo), float(hi)),
                    score_flips(truth[tsel], pred[psel], tol, scale)))
    return out


def burst_report(
    truth: np.ndarray,
    pred: np.ndarray,
    bursts,
    pad: float = 0.0,
    tol: float = DEFAULT_TOL,
    scale: float = DEFAULT_SCALE,
) -> list[dict]:
    """Per-burst efficiency and timing bias against the tunnel count.

    ``bursts`` is the list of :class:`~qpd.simulator.BurstTruth` records from
    the simulation. Each entry of the result carries ``n_qp`` (the number of
    quasiparticle tunnels the burst injected), ``n_truth`` (how many flips
    actually fall in the window, which is what limits what is recoverable), and
    the windowed :class:`FlipScore`.
    """
    rows = []
    for b in bursts:
        lo, hi = b.t_start - pad, b.t_end + pad
        tsel = (truth >= lo) & (truth < hi)
        psel = (pred >= lo - tol) & (pred < hi + tol)
        s = score_flips(truth[tsel], pred[psel], tol, scale)
        span = max(hi - lo, 1e-12)
        rows.append({
            "t_arrival": float(b.t_arrival), "n_qp": int(b.n_qp),
            "t_start": float(lo), "t_end": float(hi), "span_s": float(span),
            "flip_rate_hz": s.n_truth / span, "score": s,
        })
    return rows
