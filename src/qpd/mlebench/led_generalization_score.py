"""Frozen generalization scorer C1; independent of historical scorer B.

Exact optional one-to-one matching via sparse ordered maximum-saving DP.
Convex dead-zone quadratic costs admit an uncrossed optimum. Edges at the
cutoff have zero savings and are deterministically left unmatched.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import led_window_excess as window

VERSION = 'generalization_C1'
FP_COST, FN_COST = 3., 1.
COMPONENT_WEIGHTS = {'broad': .6, 'nominal': .25, 'stress': .15}
MAX_CANDIDATE_EDGES = 5_000_000


def timing_settings(sample_rate_hz):
    if not np.isfinite(sample_rate_hz) or sample_rate_hz <= 0:
        raise ValueError('sample_rate_hz must be finite and positive')
    tolerance = max(50e-6, 2. / sample_rate_hz)
    return tolerance, max(.002, 4*tolerance)


def _times(values):
    result = np.asarray(values, dtype=float)
    if result.ndim != 1 or not np.isfinite(result).all():
        raise ValueError('Event times must be finite one-dimensional arrays')
    if np.any(np.diff(result) < 0):
        raise ValueError('Sorted event times required')
    return result


def match_events(truth, predicted, sample_rate_hz):
    """Return original sorted-array indices and signed residuals.

    O(E log(n_pred)) time and O(E+n_pred) storage, E eligible edges, never a
    dense truth-by-prediction matrix. Excessive inputs fail explicitly rather
    than approximating the metric. Exact primary ties retain the first state.
    """
    t, p = _times(truth), _times(predicted)
    tolerance, cutoff = timing_settings(sample_rate_hz)
    left = np.searchsorted(p, t-cutoff, side='right')
    right = np.searchsorted(p, t+cutoff, side='left')
    edges = int(np.sum(right-left))
    if edges > MAX_CANDIDATE_EDGES:
        raise ValueError(f'Scorer resource limit: {edges} candidate edges; no approximate score returned')
    tree_score = np.zeros(len(p)+1)
    tree_node = np.full(len(p)+1, -1, dtype=np.int64)
    parents = np.empty(edges, dtype=np.int64)
    truth_index = np.empty(edges, dtype=np.int64)
    pred_index = np.empty(edges, dtype=np.int64)
    node_scores = np.empty(edges)
    node = 0
    for i, (lo, hi) in enumerate(zip(left, right)):
        start = node
        # Query all edges in one truth row before updates: cannot reuse truth.
        for j in range(int(lo), int(hi)):
            best, parent, k = 0., -1, j
            while k:
                if tree_score[k] > best:
                    best, parent = tree_score[k], tree_node[k]
                k -= k & -k
            distance = abs(p[j]-t[i])
            saving = 4. - 4.*(max(distance-tolerance, 0.)/(cutoff-tolerance))**2
            parents[node], truth_index[node], pred_index[node] = parent, i, j
            node_scores[node] = best + max(saving, 0.)
            node += 1
        for edge in range(start, node):
            k = int(pred_index[edge])+1
            while k <= len(p):
                if node_scores[edge] > tree_score[k]:
                    tree_score[k], tree_node[k] = node_scores[edge], edge
                k += k & -k
    best, last, k = 0., -1, len(p)
    while k:
        if tree_score[k] > best:
            best, last = tree_score[k], tree_node[k]
        k -= k & -k
    ti, pi = [], []
    while last >= 0:
        ti.append(truth_index[last]); pi.append(pred_index[last])
        last = parents[last]
    ti, pi = np.asarray(ti[::-1], dtype=int), np.asarray(pi[::-1], dtype=int)
    return ti, pi, p[pi]-t[ti]


def score_record(truth, predicted, sample_rate_hz, duration_s=None):
    scores, _ = _score_and_matches(truth, predicted, sample_rate_hz, duration_s)
    return scores


def _score_and_matches(truth, predicted, sample_rate_hz, duration_s):
    t, p = _times(truth), _times(predicted)
    if duration_s is not None:
        if not np.isfinite(duration_s) or duration_s <= 0:
            raise ValueError('duration_s must be finite and positive')
        if any(np.any((v < 0) | (v >= duration_s)) for v in (t,p)):
            raise ValueError('Events must be in [0,duration_s)')
    tolerance, cutoff = timing_settings(sample_rate_hz)
    ti, pi, residual = match_events(t,p,sample_rate_hz)
    d = np.abs(residual)
    timing = float(np.sum(4.*(np.maximum(d-tolerance,0)/(cutoff-tolerance))**2))
    precise = int(np.sum(d <= tolerance))
    fp, fn = len(p)-len(pi), len(t)-len(ti)
    scores = dict(n_truth=len(t), n_pred=len(p), n_matches=len(ti),
        n_precise=precise, n_fixed_50us=int(np.sum(d <= 50e-6)),
        n_late=len(ti)-precise, false_positives=fp, false_negatives=fn,
        timing_cost=timing, total_cost=3*fp+fn+timing,
        duration_s=duration_s, precise_s=tolerance, cutoff_s=cutoff)
    return scores, pi


def aggregate_records(records):
    records = list(records)
    keys = ('n_truth','n_pred','n_matches','n_precise','n_fixed_50us','n_late',
            'false_positives','false_negatives','timing_cost','total_cost')
    result = {key: sum(r[key] for r in records) for key in keys}
    durations = [r.get('duration_s') for r in records]
    exposure = sum(durations) if durations and all(d is not None for d in durations) else None
    result.update(n_records=len(records), duration_s=exposure,
        error=result['total_cost']/result['n_truth'] if result['n_truth'] else None,
        cost_per_second=result['total_cost']/exposure if exposure else None,
        false_positives_per_second=result['false_positives']/exposure if exposure else None)
    return result


def grade_generalization(submission, answers, manifest, *, native_submission=None):
    """DataFrames: events(record_id,t); manifest acquisition + private strata.

    Manifest requires record_id,duration_s,sample_rate_hz,component,zero_led.
    Entire missing/empty-truth strata invalidate primary score, never silently
    change weights. Individual empty recordings remain in pooled numerators.

    The supplementary LED report uses public first_pulse_time_s/on_time_s/off_time_s.
    Provide pulse_condition labels to identify design-matched active/control panels
    and a native HMM-SF event table to compare already frozen operating points.
    """
    required = {'record_id','duration_s','sample_rate_hz','component','zero_led'}
    if not isinstance(manifest,pd.DataFrame) or not manifest.columns.is_unique or not required <= set(manifest):
        raise ValueError('Invalid manifest schema')
    if manifest.empty or manifest.record_id.duplicated().any() or any(not isinstance(x,str) or not x for x in manifest.record_id):
        raise ValueError('Manifest requires unique nonempty record IDs')
    if not manifest.component.isin(COMPONENT_WEIGHTS).all() or any(not isinstance(x,(bool,np.bool_)) for x in manifest.zero_led):
        raise ValueError('Invalid component or zero_led')
    indexed = manifest.set_index('record_id')
    event_maps = []
    tables = [submission, answers]
    if native_submission is not None:
        tables.append(native_submission)
    for table in tables:
        if not isinstance(table,pd.DataFrame) or not table.columns.is_unique or set(table) != {'record_id','t'}:
            raise ValueError('Events require exactly record_id,t')
        if any(not isinstance(x,str) for x in table.record_id) or not set(table.record_id) <= set(indexed.index):
            raise ValueError('Unknown/invalid event record ID')
        if len(table) and (not pd.api.types.is_numeric_dtype(table.t) or pd.api.types.is_bool_dtype(table.t) or pd.api.types.is_complex_dtype(table.t)):
            raise ValueError('Event times must be real numeric values')
        event_maps.append({rid: np.sort(g.t.to_numpy(dtype=float)) for rid,g in table.groupby('record_id')})
    preds, truth = event_maps[:2]
    has_timing = window.TIMING_COLUMNS <= set(manifest)
    has_condition = 'pulse_condition' in manifest
    if has_condition and any(not isinstance(v, str) or not v.strip() for v in manifest.pulse_condition):
        raise ValueError('Expected nonempty pulse_condition labels')
    records = []
    pulse_records = []
    conditions = {}
    for row in manifest.itertuples(index=False):
        actual = truth.get(row.record_id, [])
        predicted = preds.get(row.record_id, [])
        r, matched = _score_and_matches(actual, predicted, row.sample_rate_hz, row.duration_s)
        r.update(record_id=row.record_id,component=row.component,zero_led=bool(row.zero_led))
        records.append(r)
        if has_timing:
            windows = window.pulse_windows(row)
            counts = dict(record_id=row.record_id, zero_led=bool(row.zero_led), pulses=len(windows[0]),
                          truth=window.window_count(actual, windows),
                          candidate=window.method_counts(predicted, matched, windows))
            if native_submission is not None:
                baseline = event_maps[2].get(row.record_id, [])
                _, baseline_matches = _score_and_matches(actual, baseline, row.sample_rate_hz, row.duration_s)
                counts['native'] = window.method_counts(baseline, baseline_matches, windows)
            pulse_records.append(counts)
            if has_condition:
                conditions.setdefault(row.pulse_condition, []).append(counts)
    strata, invalid = [], []
    for component, weight in COMPONENT_WEIGHTS.items():
        for zero_led, fraction in ((True,.1),(False,.9)):
            selected = [r for r in records if r['component']==component and r['zero_led']==zero_led]
            agg = aggregate_records(selected)
            if agg['error'] is None:
                invalid.append(f'{component}/zero_led={zero_led}: missing or zero-truth stratum')
            strata.append(dict(component=component,zero_led=zero_led,weight=weight*fraction,**agg))
    pulse_report = window.summarize_conditions(conditions, native_submission is not None)
    pulse_report.update(by_record=pulse_records, method_count_fields=['detections','false_positives'])
    if not has_timing:
        pulse_report['unavailable_reason'] = 'Public pulse timing missing'
    elif not has_condition:
        pulse_report['unavailable_reason'] = 'Matched active/control pulse_condition labels missing'
    return dict(version=VERSION,primary_valid=not invalid,invalid_reason=invalid,led_window_3ms=pulse_report,
        error=sum(s['weight']*s['error'] for s in strata) if not invalid else None,
        by_stratum=strata,by_record=records,pooled=aggregate_records(records),
        empty_truth=aggregate_records([r for r in records if not r['n_truth']]),
        settings=dict(fp_cost=3.,fn_cost=1.,precise='max(50us,2/fs)',cutoff='max(2ms,4*precise)',
                      timing='4*((abs(dt)-precise)+/(cutoff-precise))^2',aggregation='fixed six-stratum truth-normalized'))
