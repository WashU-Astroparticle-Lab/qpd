"""Companion to C1: LED-window count sensitivity at native false-positive budgets."""
import numpy as np

WINDOW_S = 0.003
TIMING_COLUMNS = {"first_pulse_time_s", "on_time_s", "off_time_s"}


def pulse_windows(row):
    first, on, off = row.first_pulse_time_s, row.on_time_s, row.off_time_s
    if not np.isfinite([first, on, off]).all() or first < 0 or on <= 0 or off <= 0:
        raise ValueError("Invalid public pulse timing")
    period = on + off
    if period < WINDOW_S:
        raise ValueError("Overlapping 3-ms pulse windows")
    # Complete half-open windows only. Include a final window ending at duration.
    starts = first + period * np.arange(max(0, int(np.ceil((row.duration_s - first) / period))))
    ends = starts + WINDOW_S
    # A pulse ending on the acquisition boundary can round one ULP above it.
    complete = ends <= np.nextafter(row.duration_s, np.inf)
    return starts[complete], ends[complete]


def window_count(sorted_times, windows):
    starts, ends = windows
    return int((np.searchsorted(sorted_times, ends) - np.searchsorted(sorted_times, starts)).sum())


def method_counts(predicted, matched_indices, windows):
    detected = window_count(predicted, windows)
    matched = window_count(np.asarray(predicted)[matched_indices], windows)
    return detected, detected - matched


def arm_counts(records, method):
    exposed = [r for r in records if r["pulses"]]
    pulses = sum(r["pulses"] for r in exposed)
    truth = sum(r["truth"] for r in exposed)
    detections = sum(r[method][0] for r in exposed)
    false_positives = sum(r[method][1] for r in exposed)
    mean = detections / pulses if pulses else None
    n = len(exposed)
    variance = None
    if n >= 2:
        # Variance of a ratio of totals, clustered by independent record.
        residuals = [pulses*r[method][0] - detections*r["pulses"] for r in exposed]
        variance = n / (n - 1) * sum(x*x for x in residuals) / pulses**4
    return dict(records=len(records), exposed_records=n, pulses=pulses, truth=truth,
                detections=detections, false_positives=false_positives, detections_per_pulse=mean,
                truth_per_pulse=truth/pulses if pulses else None,
                fp_per_100_pulses=100*false_positives/pulses if pulses else None,
                variance_of_mean=variance)


def method_summary(records, method):
    active = arm_counts([r for r in records if not r["zero_led"]], method)
    zero = arm_counts([r for r in records if r["zero_led"]], method)
    excess = truth = se = snr = None
    if active["pulses"] and zero["pulses"]:
        excess = active["detections_per_pulse"] - zero["detections_per_pulse"]
        truth = active["truth_per_pulse"] - zero["truth_per_pulse"]
    if active["variance_of_mean"] is not None and zero["variance_of_mean"] is not None:
        se = float(np.sqrt(active["variance_of_mean"] + zero["variance_of_mean"]))
        if se > 0:
            snr = excess / se
    return dict(active=active, zero_led=zero, excess_per_pulse=excess,
                true_excess_per_pulse=truth, excess_bias_per_pulse=None if truth is None else excess-truth,
                record_cluster_se=se, snr=snr)


def summarize_conditions(groups, has_native):
    results = {}
    for name, records in groups.items():
        candidate = method_summary(records, "candidate")
        native = method_summary(records, "native") if has_native else None
        passes = {}
        for arm in ("active", "zero_led"):
            passes[arm] = (candidate[arm]["false_positives"] <= native[arm]["false_positives"]
                           if native is not None and candidate[arm]["pulses"] else None)
        within_budget = all(passes.values()) if None not in passes.values() else None
        results[name] = dict(candidate=candidate, native=native, passes_fp_budget=passes,
                             both_budgets_pass=within_budget,
                             snr_at_native_fp_budget=candidate["snr"] if within_budget else None)
    return dict(version="led_excess_3ms_v1", window_s=WINDOW_S, by_condition=results)
