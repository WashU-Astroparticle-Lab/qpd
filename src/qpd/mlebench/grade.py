"""Grader for the MLE-bench QPD datasets.

Both submissions and answer keys are long-format event tables with columns
``[chunk_id, event_type, t, value]`` (see :mod:`qpd.mlebench.generate`).

Grading is per chunk and per event type:

1. **Matching** -- predicted events are matched to truth events *optimally*
   (minimum total timing error, Hungarian / ``linear_sum_assignment``) subject
   to a per-type time tolerance. To stay tractable on dense flip trains
   (~thousands per chunk) the tolerance-gated bipartite graph is split into
   connected components and the assignment is solved per component; this is
   exact because no edge crosses a component boundary.

2. **Continuous timing score** -- each matched event earns partial credit
   ``q = exp(-(dt / timing_scale)^2)`` from its timing residual ``dt``. A
   *soft* F1 folds this into detection: ``soft_TP = sum(q)`` over matches,
   ``soft_precision = soft_TP / n_pred``, ``soft_recall = soft_TP / n_truth``.
   A perfectly-timed, complete match scores 1; false positives/negatives and
   loose timing both pull it down.

The scalar :func:`grade` returns the mean soft-F1 across the event types that
occur (in either the key or the submission), lightly penalised by the
normalised value error (fractional charge / QP count) of the matched events.
:func:`grade_report` returns the full per-type breakdown and, optionally, the
per-event match records.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

# Default time-matching tolerances [s], per event type (hard gate).
DEFAULT_TOLERANCES = {
    "flip": 0.5e-3,
    "charge": 5.0e-3,
    "burst": 5.0e-3,
}
# Timing kernel scale [s] per type. Default = tolerance / 2 (a match right at
# the gate edge then earns exp(-4) ~ 0.018; a tight match ~ 1).
_DEFAULT_SCALE_FRACTION = 0.5
# Scale for normalising value errors into the [0, 1] penalty.
_VALUE_SCALE = {"charge": 0.5, "burst": 15.0}
_EVENT_TYPES = ("flip", "charge", "burst")


@dataclass
class MatchResult:
    """Optimal one-to-one matching of a truth and a prediction time set."""

    n_truth: int
    n_pred: int
    matched_truth_idx: np.ndarray = field(default_factory=lambda: np.empty(0, int))
    matched_pred_idx: np.ndarray = field(default_factory=lambda: np.empty(0, int))
    # Signed timing residuals dt = pred - truth for each match [s].
    dt: np.ndarray = field(default_factory=lambda: np.empty(0, float))

    @property
    def tp(self) -> int:
        return int(self.matched_truth_idx.size)

    @property
    def fp(self) -> int:
        return self.n_pred - self.tp

    @property
    def fn(self) -> int:
        return self.n_truth - self.tp

    @property
    def f1(self) -> float:
        tp = self.tp
        denom = 2 * tp + self.fp + self.fn
        if denom == 0:
            return 1.0
        return 2 * tp / denom


def _components(nt: int, npred: int, edges: list[tuple[int, int]]):
    """Union-find over truth nodes [0,nt) and pred nodes [0,npred); group by root.

    Returns a list of (truth_idx_list, pred_idx_list); only nodes that appear in
    at least one candidate edge are included (others are necessarily unmatched).
    """
    parent = list(range(nt + npred))

    def find(x):
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:
            parent[x], x = root, parent[x]
        return root

    for i, j in edges:
        ri, rj = find(i), find(nt + j)
        if ri != rj:
            parent[ri] = rj

    groups: dict[int, tuple[list, list]] = {}
    for i, j in edges:
        r = find(i)
        if r not in groups:
            groups[r] = ([], [])
        # Each edge endpoint is added; dedup afterwards.
        groups[r][0].append(i)
        groups[r][1].append(j)
    return [
        (sorted(set(ti)), sorted(set(pj))) for ti, pj in groups.values()
    ]


def match_times(
    truth_t: np.ndarray, pred_t: np.ndarray, tol: float
) -> MatchResult:
    """Optimal (min total |dt|) one-to-one matching within ``tol``.

    The tolerance-gated bipartite graph is decomposed into connected components
    and ``linear_sum_assignment`` is solved per component. Exact and fast: each
    component is a local cluster of events, so the assignments are tiny.
    """
    truth_t = np.sort(np.asarray(truth_t, dtype=float))
    pred_t = np.sort(np.asarray(pred_t, dtype=float))
    nt, npred = truth_t.size, pred_t.size
    if nt == 0 or npred == 0:
        return MatchResult(nt, npred)

    # Candidate edges (|dt| <= tol) via a sliding window over sorted arrays.
    edges: list[tuple[int, int]] = []
    j0 = 0
    for i, tt in enumerate(truth_t):
        while j0 < npred and pred_t[j0] < tt - tol:
            j0 += 1
        j = j0
        while j < npred and pred_t[j] <= tt + tol:
            edges.append((i, j))
            j += 1
    if not edges:
        return MatchResult(nt, npred)

    mt: list[int] = []
    mp: list[int] = []
    mdt: list[float] = []
    big = tol * 10.0 + 1.0
    for ti, pj in _components(nt, npred, edges):
        if len(ti) == 1 and len(pj) == 1:
            d = pred_t[pj[0]] - truth_t[ti[0]]
            if abs(d) <= tol:
                mt.append(ti[0])
                mp.append(pj[0])
                mdt.append(d)
            continue
        cost = np.abs(
            truth_t[np.asarray(ti)][:, None] - pred_t[np.asarray(pj)][None, :]
        )
        gated = np.where(cost <= tol, cost, big)
        rows, cols = linear_sum_assignment(gated)
        for a, b in zip(rows, cols):
            if cost[a, b] <= tol:
                mt.append(ti[a])
                mp.append(pj[b])
                mdt.append(pred_t[pj[b]] - truth_t[ti[a]])

    order = np.argsort(mt) if mt else np.empty(0, int)
    return MatchResult(
        n_truth=nt,
        n_pred=npred,
        matched_truth_idx=np.asarray(mt, int)[order],
        matched_pred_idx=np.asarray(mp, int)[order],
        dt=np.asarray(mdt, float)[order],
    )


def _circular_abs(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Absolute difference of fractional charges on a period-1 circle."""
    d = (np.asarray(a) - np.asarray(b) + 0.5) % 1.0 - 0.5
    return np.abs(d)


def _subset(df: pd.DataFrame, chunk_id, event_type) -> pd.DataFrame:
    return df[(df["chunk_id"] == chunk_id) & (df["event_type"] == event_type)]


def grade_report(
    submission: pd.DataFrame,
    answers: pd.DataFrame,
    tolerances: dict | None = None,
    timing_scales: dict | None = None,
    return_matches: bool = False,
) -> dict:
    """Full grading breakdown. See module docstring.

    Parameters
    ----------
    return_matches : bool
        If True, include a ``matches`` DataFrame with one row per matched event
        ``[chunk_id, event_type, t_truth, t_pred, dt, quality]`` for diagnostics.
    """
    tol = {**DEFAULT_TOLERANCES, **(tolerances or {})}
    scale = {k: _DEFAULT_SCALE_FRACTION * v for k, v in tol.items()}
    scale.update(timing_scales or {})

    for col in ("chunk_id", "event_type", "t"):
        if col not in answers.columns:
            raise ValueError(f"answers missing column {col!r}")
    if submission is None or len(submission) == 0:
        submission = pd.DataFrame(columns=["chunk_id", "event_type", "t", "value"])

    chunk_ids = sorted(answers["chunk_id"].unique())
    present = set(answers["event_type"].unique()) | set(
        submission["event_type"].unique()
    )
    event_types = [et for et in _EVENT_TYPES if et in present]

    match_records: list[dict] = []
    per_type: dict[str, dict] = {}
    for et in event_types:
        s = scale[et]
        n_truth_tot = n_pred_tot = 0
        hard_tp = 0
        soft_tp = 0.0
        timing_sq: list[float] = []
        value_err: list[float] = []
        for cid in chunk_ids:
            tr = _subset(answers, cid, et).sort_values("t")
            pr = _subset(submission, cid, et).sort_values("t")
            tr_t = tr["t"].to_numpy()
            pr_t = pr["t"].to_numpy()
            n_truth_tot += tr_t.size
            n_pred_tot += pr_t.size
            m = match_times(tr_t, pr_t, tol[et])
            if m.tp == 0:
                continue
            hard_tp += m.tp
            q = np.exp(-((m.dt / s) ** 2)) if s > 0 else (m.dt == 0).astype(float)
            soft_tp += float(q.sum())
            timing_sq.extend((m.dt ** 2).tolist())
            if et in ("charge", "burst"):
                tv = tr["value"].to_numpy()[m.matched_truth_idx]
                pv = pr["value"].to_numpy()[m.matched_pred_idx]
                if et == "charge":
                    value_err.extend(_circular_abs(pv, tv).tolist())
                else:
                    value_err.extend(np.abs(pv - tv).tolist())
            if return_matches:
                ttt = tr_t[m.matched_truth_idx]
                ptt = pr_t[m.matched_pred_idx]
                for a in range(m.tp):
                    match_records.append(
                        {"chunk_id": cid, "event_type": et,
                         "t_truth": float(ttt[a]), "t_pred": float(ptt[a]),
                         "dt": float(m.dt[a]), "quality": float(q[a])}
                    )

        # Hard detection (diagnostic).
        hard_prec = hard_tp / n_pred_tot if n_pred_tot else (1.0 if n_truth_tot == 0 else 0.0)
        hard_rec = hard_tp / n_truth_tot if n_truth_tot else 1.0
        hard_f1 = (2 * hard_prec * hard_rec / (hard_prec + hard_rec)
                   if (hard_prec + hard_rec) else 0.0)
        # Continuous timing-aware soft F1 -- THE per-type score.
        soft_prec = soft_tp / n_pred_tot if n_pred_tot else (1.0 if n_truth_tot == 0 else 0.0)
        soft_rec = soft_tp / n_truth_tot if n_truth_tot else 1.0
        timing_f1 = (2 * soft_prec * soft_rec / (soft_prec + soft_rec)
                     if (soft_prec + soft_rec) else
                     (1.0 if n_truth_tot == 0 and n_pred_tot == 0 else 0.0))

        per_type[et] = {
            "n_truth": n_truth_tot,
            "n_pred": n_pred_tot,
            "hard_tp": hard_tp,
            "hard_precision": hard_prec,
            "hard_recall": hard_rec,
            "hard_f1": hard_f1,
            "soft_tp": soft_tp,
            "soft_precision": soft_prec,
            "soft_recall": soft_rec,
            "timing_f1": timing_f1,
            "mean_match_quality": (soft_tp / hard_tp) if hard_tp else None,
            "timing_rmse_s": float(np.sqrt(np.mean(timing_sq))) if timing_sq else None,
            "value_mae": float(np.mean(value_err)) if value_err else None,
            "tolerance_s": tol[et],
            "timing_scale_s": s,
        }

    # Overall: mean timing_f1 across present types, penalised by value error.
    if per_type:
        mean_f1 = float(np.mean([v["timing_f1"] for v in per_type.values()]))
        penalties = [
            min(1.0, v["value_mae"] / _VALUE_SCALE[et])
            for et, v in per_type.items()
            if v["value_mae"] is not None and et in _VALUE_SCALE
        ]
        value_penalty = float(np.mean(penalties)) if penalties else 0.0
        overall = mean_f1 * (1.0 - 0.5 * value_penalty)
    else:
        mean_f1 = overall = value_penalty = 0.0

    report = {
        "overall_score": overall,
        "mean_timing_f1": mean_f1,
        "value_penalty": value_penalty,
        "n_chunks": len(chunk_ids),
        "event_types": event_types,
        "per_type": per_type,
        "tolerances": tol,
        "timing_scales": scale,
    }
    if return_matches:
        report["matches"] = pd.DataFrame(
            match_records,
            columns=["chunk_id", "event_type", "t_truth", "t_pred", "dt", "quality"],
        )
    return report


def grade(
    submission: pd.DataFrame,
    answers: pd.DataFrame,
    tolerances: dict | None = None,
    timing_scales: dict | None = None,
) -> float:
    """MLE-bench-style scalar grade in [0, 1] (higher is better)."""
    return float(
        grade_report(submission, answers, tolerances, timing_scales)[
            "overall_score"
        ]
    )
